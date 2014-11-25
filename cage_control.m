classdef cage_control < hgsetget
    %CAGE_CONTROL class to control the Helmholtz cage
    %   this class provides an interface to the Helmholtz cage power
    %   supplies and magnetometer
    %
    %cage_control methods:
    %   BtoI            - convert magnetic field vector to current vector
    %	ItoB            - convert current vector to magnetic field vector
    %	calibrate       - generate a new calibration matrix
    %	cal_test        - test calibration matrix
    %   loadCal         - load calibration matrix
    %   saveCal         - save calibration matrix
    %
    %cage_control properties:
    %   T       - Calibration matrix, used to convert from magnetic field to power supply current
    %   Is      - Current set point
    %   Bs      - Magnetic Field Set Point (calculated from Is)
    %	Bm      - Magnetic Field Reading

    

    
    properties
        T=[2 0 0 0;0 2 0 0;0 0 2 0;0 0 0 1];    %Calibration matrix, used to convert from magnetic field to power supply current
        Is=[0 0 0];                             %Current set point
    end
    
    properties(Dependent = true)
        Bs                                      %Magnetic Field Set Point (calculated from Is)
    end
    
    properties(SetAccess=private,Dependent = true,Hidden)
        Bm                                      %Magnetic Field Reading
    end
    
    properties(GetAccess=private,SetAccess=private,Hidden)
        x_obj                                   
        y_obj
        z_obj
        mag
        dir
    end
    
    properties(SetAccess = private,GetAccess = private,Dependent = true,Hidden)
        psObj
        allObj
    end
    
    properties(Constant,Hidden)
            %power supply addresses
            x_addr=[1 2];   %addresses for x axis
            y_addr=[3 4];   %addresses for y axis
            z_addr=[5 6];   %addresses for z axis
            addr_names={'X1','X2','Y1','Y2','Z1','Z2'};
    end

    %static methods, called by refernecing classname rather than class
    %instance
    methods(Static)
        function session=open_mag()
            % ai=open_mag()
            % used to open an analog input interface to the magnetometer

            
            %open data acquisition card
            session = daq.createSession('ni');
            %names for mag channels
            names={'X-mag','Y-mag','Z-mag'};
            %add magnetometer X, Y, and Z channels
            for k=0:2
                c=session.addAnalogInputChannel('Helmholtz',['ai' num2str(k)],'Voltage');
                set(c,'Coupling','DC');
                set(c,'TerminalConfig','SingleEnded');
                set(c,'Range',[-5 5]);
                set(c,'Name',names{k+1});
            end
           
            %set sample rate to 10 sps
            set(session,'Rate',10) 
            %set scan time to 1 sec
            set(session,'DurationInSeconds',1);
        end
        function close_mag(session)
            %close_mag(session)
            %close the magnetomitor interface and stop all data acquisition
            
            %stop and close data acquisition
            session.stop();
            delete(session);    
        end
    end
    
    methods
        function obj=cage_control(mag)
            if(nargin==0)
                mag=true;
            end

            if(is_gpib_open([obj.x_addr obj.y_addr obj.z_addr]))
                resp=questdlg('One or more power supplies are currently open. Close to proceed?','Helmholtz Cage Init','Yes','No','Yes');
                if isequal(resp,'Yes')
                    %close all used supplies
                    force_gpib_close([obj.x_addr obj.y_addr obj.z_addr]);
                else
                    %else throw an error
                    error('HELMHOLTZ:invalidCurr','Unable to open power supplies.');
                end
            end

            %open GPIB interface to X-axis power supplies
            obj.x_obj=open_supplys(obj.x_addr,'X');
            %open GPIB interface to Y-axis power supplies
            obj.y_obj=open_supplys(obj.y_addr,'Y');
            %open GPIB interface to Z-axis power supplies
            obj.z_obj=open_supplys(obj.z_addr,'Z');

            if(mag)
                %open magnetometer
                obj.mag=obj.open_mag();
                Bm=get(obj,'Bm');
            else
                obj.mag=[];
            end

            %open digital lines for direction control
            %obj.dir=digitalio('nidaq','Dev1');
            %add lines for X,Y, and Z direction
            %addline(obj.dir,0:2,1,'out');
            %set all directions to zero
            %putvalue(obj.dir,0);
        end
        function delete(obj)
           %stop and close the power supplies
            for k=(1:length(obj.allObj))
                close_supply(obj.allObj(k));
            end
            if(isempty(obj.mag)==false)
                %close the magnetometer
                cage_control.close_mag(obj.mag);
            end
            %clear all outputs
            %putvalue(obj.dir,0)
            %delete digitial IO object
            %delete(obj.dir);
        end
        function val=get.psObj(obj)
            val={obj.x_obj,obj.y_obj,obj.z_obj};
        end
        function val=get.allObj(obj)
            val=[obj.x_obj,obj.y_obj,obj.z_obj];
        end
        function obj=set.Is(obj,I)
            if ~isvector(I) || length(I)~=3
                error('HELMHOLTZ:invalidCurr','Supply Current Vector ''Is'' must be a 3 element vector');
            end
            %make I a row vector
            I=reshape(I,1,3);
            %round to power supply resolution
            I=round(I*1e3)*1e-3;
            %set each axis separately
            for k=1:3
                %only set axis that have changed
                if I(k)~=obj.Is(k)
                    %get power supply objects
                    g=obj.psObj{k};
                    %current value
                    It=I(k);
                    neg=It<0;
                    %switch current direction for negative currents
                    %putvalue(obj.dir.Line(k),neg);
                    %get magnitude of current
                    It=abs(It);
                    %set the maximum current to 8A and give a warning
                    if  It>7
                        It=7;
                        %also change I as well
                        I(k)=7-(14*neg);
                        warning('HELMHOLTZ:invalidCurr','Supply current too large');
                    end
                    %set the currents for the power supplies
                    for kk=1:length(g)
                        fprintf(g(kk),'SOUR:CURR %f',It)
                    end
                end
            end
            obj.Is=I;
        end
        function I=BtoI(obj,B)
            %I=BtoI(B)
            %returns the reacquired current vector to produce the given
            %magnetic field vector 'B'.
            
            s=size(B);
            if(length(s)~=2) 
                error('Argument must be a vector')
            end
            if(s(1)==1)
                B=B';
                s=size(B);
            end

            if(s(2)~=1)
                error('Argument must be a vector')
            end

            len=length(B);
            if(len<4)
                B=[B;ones(4-len,1)];
            end

            %do the actual transformation
            I=obj.T*B;
            %extract current
            I=I(1:3)';
        end
        function B=ItoB(obj,I)
            %B=ItoB(I)
            %returns a magnetic field vector for a given current vector
            
            s=size(I);
            if(length(s)~=2) 
                error('Argument must be a vector')
            end
            if(s(1)==1)
                I=I';
                s=size(I);
            end

            if(s(2)~=1)
                error('Argument must be a vector')
            end

            len=length(I);
            if(len<4)
                I=[I;ones(4-len,1)];
            end

            %do the actual transformation
            B=obj.T^-1*I;
            %extract current
            B=B(1:3)';
        end
        function obj=set.Bs(obj,Bs)
            %set the magnetic field of the Helmholtz cage to Bs
            
            %check to make sure that the magnetic field is a vector before
            %proceeding
            if ~isvector(Bs) || length(Bs)~=3
                error('HELMHOLTZ:invalidField','Magnetic Field Vector ''Bs'' must be a 3 element vector');
            end
            %set Helmholtz cage field to Bs
            I=obj.BtoI(Bs);
            %set currents to calculated values
            obj.Is=I;
        end
        function test(obj)
            %check to see if connections to coils are working by running
            %test currents through the coils
            
            %save old Is
            I_save=obj.Is;
            %set all currents to 5A
            obj.Is=5*ones(1,3);
            %initialize arrays
            V=zeros(size(obj.allObj));
            I=zeros(size(obj.allObj));
            %test each supply
            for k=(1:length(obj.allObj))
                %measure voltage
                fprintf(obj.allObj(k),'MEAS:VOLT?');
                %get result
                V(k)=str2double(fgetl(obj.allObj(k)));
                %measure current
                fprintf(obj.allObj(k),'MEAS:CURR?');
                %get result
                I(k)=str2double(fgetl(obj.allObj(k)));
            end
            %calculate resistance
            R=V./I;
            %check if resistance is greater than 20 ohms
            open_coil=R>20;
            %check for any open coils
            if(any(open_coil))
                %restore Is
                obj.Is=I_save;
                %get number of open coils
                num=sum(open_coil);
                %check if number is one
                if(num==1)
                    %give error stating which coil is not working
                    error('HELMHOLTZ:badConnection','The %s coil has a bad connection',obj.addr_names{open_coil});
                else
                    %give an error stating how many coils are not working
                    error('HELMHOLTZ:badConnection','There are %i coils with bad connections',num);
                end
            end
            %restore Is
            obj.Is=I_save;
        end
        function err=calibrate(obj)
            % err=calibrate()
            % calibrate() will generate the calibration matrix for a
            % cage_control object by setting the Helmholtz cage power
            % supply currents and reading the field with the magnetometer
            % calibrate() returns the error given by the least squares fit
            % of the data. if the calibration was aborted by the user
            % calibrate() returns inf.
            
            %make shure that power supplies are on
            obj.on();
            
            hP=waitbar(0,'Calibrating Helmholtz Cage Please Wait...','WindowStyle','modal',...
                'CreateCancelBtn',@prog_close,'name','Calibrating');

            %calibration current values
            %cal_cur=(0:0.1:6);
            %cal_cur=(0:0.1:3);
            %cal_cur=(0:0.5:3);
            %cal_cur=(0:0.5:3);
            cal_cur=(-1.5:1:1.5);
            len=length(cal_cur);

            Bc=zeros(len^3,3);
            Ic=zeros(len^3,3);
            n=1;
            %range of indices to iterate over for Y and Z
            %these get switched so that the maximum step is lower
            Yrng=(1:length(cal_cur));
            Zrng=(1:length(cal_cur));
            %get start time to calculate elapsed time
            tst=tic();
            %drive all combinations of currents on all axis
            for Ixi=(1:length(cal_cur))
                %set x-axis current
                obj.Is(1)=cal_cur(Ixi);
                for Iyi=Yrng
                    %set y-axis current
                    obj.Is(2)=cal_cur(Iyi);
                    for Izi=Zrng
                        %set z-axis current
                        obj.Is(3)=cal_cur(Izi);
                        %stop any previous sampling
                        stop(obj.mag)
                        %pause to let the supply settle
                        pause(0.01);
                        %enter values in array
                        Bc(n,:)=obj.Bm;
                        Ic(n,:)=[cal_cur(Ixi),cal_cur(Iyi),cal_cur(Izi)];
                        %increment array index
                        n=n+1;
                        %check progress bar handle to see if it has been canceled
                        if not(ishandle(hP))
                            break;
                        end
                        %done fraction
                        df=n/len^3;
                        %update progress bar
                        waitbar(df,hP);
                        %get time elapsed
                        Te=toc(tst);
                        %get remaining time
                        Tr=Te*(df^-1-1);
                        %update window title
                        set(hP,'name',sprintf('Calibrating (%s remaining)',sec2str(Tr)));
                    end
                    Zrng=Zrng(end:-1:1);
                    if not(ishandle(hP))
                        break;
                    end
                end
                Yrng=Yrng(end:-1:1);
                if not(ishandle(hP))
                    break;
                end
            end

            if not(ishandle(hP))
                %set err to inf to indicate failure
                err=inf;
                return;
            end

            %close progress bar
            close(hP);

            %calculate calibration transform matrix
            %setup coefficient matrix
            A=[Bc,ones(length(Bc),1)];
            %use least squares to get solution matrix
            sol=(A'*A)\(A');
            %use the least squared solution with each I vector to get the
            %transform matrix
            xhx=sol*Ic(:,1);
            xhy=sol*Ic(:,2);
            xhz=sol*Ic(:,3);
            obj.T=[(xhx)';(xhy)';(xhz)';0 0 0 1];

            ex=sqrt(mean((Ic(:,1)-A*xhx).^2));
            ey=sqrt(mean((Ic(:,2)-A*xhy).^2));
            ez=sqrt(mean((Ic(:,3)-A*xhz).^2));

            %fprintf(1,'Calibration Error X = %f\n',ex);
            %fprintf(1,'Calibration Error Y = %f\n',ey);
            %fprintf(1,'Calibration Error Z = %f\n',ez);
            err=mean([ex ey ez]);

            if(err>0.5)
                choice=questdlg('Calibration Error Large. Use Calibration anyway?','Calibration Error','Yes','No','Try Again','No');
                %handle user choice
                switch  choice
                    case 'Try Again'
                        %run calibration again
                        err=obj.calibrate();
                    case 'No'
                        %don't change calibration
                        err=inf;
                   %Yes falls through
                end
            end
    
        end
        function [max_err,mean_err]=cal_test(obj)
            % [max_err,mean_err]=cal_test()
            % test the calibration matrix of the cage_control object
            % cal_test will compare the measured vs. expected values
            % for a number of different magnetic field values then plot a
            % graph and return the mean and max error  
            
            %make shure that power supplies are on
            obj.on();
            hP=waitbar(0,'Testing Helmholtz Cage Calibration Please Wait...','WindowStyle','modal',...
                'CreateCancelBtn',@prog_close,'name','Testing');

            %set number of samples per trigger
            set(obj.mag,'SamplesPerTrigger',10);

            %calibration field values
            %xrng=(0:0.05:0.5);
            %yrng=(0:-0.05:-0.5);
            %zrng=(0:0.05:0.5);
            rng=(-0.5:0.1:0.5);
            xrng=rng;
            yrng=rng;
            zrng=rng;
            
            %preallocate for speed like matlab tells you to
            errX=zeros(length(xrng),length(yrng),length(zrng));
            errY=zeros(length(xrng),length(yrng),length(zrng));
            errZ=zeros(length(xrng),length(yrng),length(zrng));
            %maximum value for progress indicator
            maxP=length(xrng)*length(yrng)*length(zrng);

            Yrng=(1:length(yrng));
            Zrng=(1:length(zrng));

            %get start time to calculate elapsed time
            tst=tic();
            %keep track of iterations for progress
            n=1;
            %drive all combinations of field and measure error
            for Bxi=(1:length(xrng))
                %set the x-axis field
                obj.Bs(1)=xrng(Bxi); 
                for Byi=Yrng
                    %set the y-axis field
                    obj.Bs(2)=yrng(Byi); 
                    for Bzi=Zrng
                        %set the z-axis field
                        obj.Bs(3)=zrng(Bzi);     
                        %stop measurments
                        stop(obj.mag);
                        %wait for field to settle
                        pause(0.3);
                        %average measurement
                        dat=obj.Bm;
                        %compute error
                        errX(Bxi,Byi,Bzi)=dat(1)-xrng(Bxi);
                        errY(Bxi,Byi,Bzi)=dat(2)-yrng(Byi);
                        errZ(Bxi,Byi,Bzi)=dat(3)-zrng(Bzi);
                        %check to see if progress bar was closed
                        if not(ishandle(hP))
                            break;
                        end
                        %increment iteration counter
                        n=n+1;
                       %done fraction
                        df=(n)/maxP;
                        %update progress bar
                        waitbar(df,hP);
                        %get time elapsed
                        Te=toc(tst);
                        %get remaining time
                        Tr=Te*(df^-1-1);
                        %update window title
                        set(hP,'name',sprintf('Testing (%s remaining)',sec2str(Tr)));
                    end
                    %check to see if progress bar was closed
                    if not(ishandle(hP))
                        max_err=NaN;
                        mean_err=NaN;
                        return;
                    end
                    %switch Z range
                    Zrng=Zrng(end:-1:1);
                end

                %switch Y range
                Yrng=Yrng(end:-1:1);
            end


            %close progress bar
            close(hP);
            %create a figure for the data
            hFig=figure();
            hAxes=axes();
            %make plot of Field Error Data
            %min and max values
            xmax=max(xrng);
            xmin=min(xrng);
            ymax=max(yrng);
            ymin=min(yrng);
            zmax=max(zrng);
            zmin=min(zrng);
            %creat a grid for the plots
            %[cx cy cz] = meshgrid(xrng,yrng,zrng);
            [cy cx cz] = meshgrid(yrng,xrng,zrng);
            % Use daspect to set the data aspect ratio of the axes
            daspect(hAxes,[1,1,1]);
            %calculate the fiedld magnitude
            err_mag = sqrt(errX.^2 + errY.^2 + errZ.^2);
            %calculate stats
            mean_err=mean(mean(mean(err_mag)));
            max_err=max(max(max(err_mag)));
            % Create slice planes for the error data in a cube shape
            hsurfaces = slice(hAxes,cy,cx,cz,err_mag,[ymin ymax],[xmin xmax],[zmin zmax]);
            %show colorbar for plot
            colorbar('WestOutside')
            % Specify interpolated face color and do not draw edges
            %set(hsurfaces,'FaceColor','interp','EdgeColor','none')

            % Use the axis command to set the axis limits equal to the range of the data.
            axis(hAxes,'tight');
            % Orient the view to azimuth = 30 and elevation = 40.
            view(hAxes,30,40);
            % Select perspective projection to provide a more realistic looking volume 
            % using camproj:
            camproj(hAxes,'orthographic'); 
            % Zoom in on the scene a little
            camzoom(hAxes,1.1);
            % Add a light source to the right of the camera and use Phong lighting to give the 
            % slice planes a smooth, three-dimensional appearance using camlight and lighting:
            hLight=light('Parent',hAxes); 
            camlight(hLight,'right');
            lighting(hAxes,'phong');
            % Increase the value of the AmbientStrength property for each slice plane to improve 
            % the visibility of the dark blue colors: 
            set(hsurfaces,'AmbientStrength',.6);
            %label the axies
            xlabel(hAxes,'x-axis Field [Gauss]');
            ylabel(hAxes,'y-axis Field [Gauss]');
            zlabel(hAxes,'z-axis Field [Gauss]');
        end
        function Bm=get.Bm(obj)
            %read measurement
            Bm=obj.mag.inputSingleScan();
        end
        function Bs=get.Bs(obj)
            Bs=obj.ItoB(obj.Is);
        end
        function loadCal(obj,name)
            %loadCal(name)
            %load the calibration matrix from the file given by name
            %if no name is given then the default calibration is loaded
            %from calibration.cal
            
            if nargin==1
                name='calibration.cal';
            end
            obj.T=load(name);
        end
        function saveCal(obj,name)
            %saveCal(name)
            %saves the current calibration matrix in the file whit the name
            %given by name
            dlmwrite(name,obj.T);
        end   
        function obj=set.T(obj,T)
            if size(T)~=[4 4]
                error('HELMHOLTZ:invalidTransform','Transformation matrix T must be a 4x4 matrix');
            end
            obj.T=T;
        end
        function on(obj)
            %on() turn Helmholtz Cage Power supplies on
            
            %turn off the power supplies
            for k=(1:length(obj.allObj))
                fprintf(obj.allObj(k),'OUTP:STAT ON');
            end
        end
        function off(obj)
            %off() turn Helmholtz Cage power supplies off
            
            %turn off the power supplies
            for k=(1:length(obj.allObj))
                fprintf(obj.allObj(k),'OUTP:STAT OFF');
            end
        end
    end
end

%the following helper functions are used within the class but are not
%avalible externally

%open a supply and set initial current and volatge values
function objs=open_supplys(addrs,axis)
    for k=(1:length(addrs))
        %open GPIB device
        objs(k)=gpib('ni', 0,addrs(k),'Name',['Helmholtz ' axis int2str(k) ' Axis']);
        fopen(objs(k));
        %set current to zero
        fprintf(objs(k),'SOUR:CURR 0');
        %set voltage limit to 10V
        fprintf(objs(k),'SOUR:VOLT 10');
        %turn on the output
        fprintf(objs(k),'OUTP:STAT ON');
    end
end
%close and turn off a supply
function close_supply(g)
    %turn output off
    fprintf(g,'OUTP:STAT OFF')
    %close object
    fclose(g)
    %delete object
    delete(g)
end

function prog_close(src,evnt)
    %close request for progress bar    

    %close progress bar
    delete(gcf)
end

function open=is_gpib_open(addr)
%check to see if supplies are open    
    for a=addr
        inst=instrfind({'Type','PrimaryAddress','Status'},{'gpib',a,'open'});
        if(size(inst))
            open=1;
            return
        end
    end
    open=0;
end

function s=sec2str(sec)
    %convert a time in seconds to a string with minutes and hours
    %for progress usage
    min=floor(sec/60);
    sec=sec-60*min;
    hour=floor(min/60);
    min=min-hour*60;
    %check for hours
    if(hour > 0)
        %1 or more hours, use full time format
        s=sprintf('%02.0f:%02.0f:%02.0f',hour,min,sec);
    else
        if min>0
            %print seconds and minutes
            s=sprintf('%02.0f:%02.0f',min,sec);
        else
            if sec <= 2
                %not much time left
                s='a few moments';
            else
                %a few seconds left, no leading zero
                s=sprintf('%.0f sec',sec);
            end
        end
    end
end
        
%force all gpib interface objects that address devices whith given
%addresses to close
function open=force_gpib_close(addr)
    for a=addr
        inst=instrfind({'Type','PrimaryAddress','Status'},{'gpib',a,'open'});
        if(size(inst))
            delete(inst);
        end
    end
end

