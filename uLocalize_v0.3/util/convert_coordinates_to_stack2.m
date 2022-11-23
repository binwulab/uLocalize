function stack = convert_coordinates_to_stack2(varargin)

%convert_coordinates_to_stack([x y z I],nx,ny,nz)
%convert_coordinates_to_stack([x y z],nx,ny,nz)
%convert_coordinates_to_stack([x y z I],nx,ny,nz,valmode)

%convert_coordinates_to_stack([x y I],nx,ny)
%convert_coordinates_to_stack([x y],nx,ny)
%convert_coordinates_to_stack([x y I],nx,ny,valmode)

%valmode optional argument
%valmode = 'value' (default): I build a stack with the values corresponding to the last
%column of the array [x y z I] or [x y I]

%valmode = 'ones': I put ones in the stack at the positions indicated by the
%input array [x y z I] or [x y I]
if nargin == 0
    stack = 0;
    return;
end

if nargin >=1
    coor = varargin{1};
    %% 3d case
    if nargin == 5 || (nargin == 4 && ~ischar(varargin{4}) )  %3d case
        nx = varargin{2}; ny = varargin{3}; nz = varargin{4};
        if nargin ==5
            valmode = varargin{5};
        else
            if size(coor,2) >= 4
                valmode = 'value';
            elseif size(coor,2) == 3
                valmode = 'ones';
            else
                disp(['3d convert_coordinates_to_stack: input array has the wrong column number: ' num2str(size(coor,2))]);
            end
        end

        stack = fill_stack_3d(coor,nx,ny,nz,valmode);
        return;

    %% 2d case
    elseif nargin == 3 || (nargin == 4 && ischar(varargin{4}) )  
        nx = varargin{2}; ny = varargin{3};
        if nargin ==4
            valmode = varargin{4};
        else

            if size(coor,2) >= 3
                valmode = 'value';
            elseif size(coor,2) == 2
                valmode = 'ones';
            else
                disp(['2d convert_coordinates_to_stack: input array has the wrong column number: ' num2str(size(coor,2))]);
            end
        end

        stack = fill_stack_2d(coor,nx,ny,valmode);
        return

    %% wrong number/format of arguments    
    else        
        return
    end
end
clear('valmode','nx','ny','nz');
end


function stack = fill_stack_3d(coor,nx,ny,nz,valmode)
% fill in a 3D stack (1:nx,1:ny,1:nz) with zeros except at positions indicated by coor = [x y z I] 
%or [x y z].
%if input positions falls outside of the size of the array, I change the coordinate
%to zero or nx/ny/nz
    
    if size(coor,2) == 3
        stack = zeros(nx,ny,nz,'uint8');
    elseif size(coor,2) == 4
        stack = zeros(nx,ny,nz,class(coor));
    end

    if ndims(coor)<2
        return
    end

    %I impose that coordinates fall within the destination array
    coor(:,1) = max(1,min(uint16(ceil(coor(:,1))),nx));
    coor(:,2) = max(1,min(uint16(ceil(coor(:,2))),ny));
    coor(:,3) = max(1,min(uint16(ceil(coor(:,3))),nz));

    for i = 1:size(coor,1)
        if size(coor,2)>=4 && strcmp(valmode,'value')            
            stack(coor(i,1),coor(i,2),coor(i,3))=coor(i,4);
        elseif strcmp(valmode,'ones')  
            stack(coor(i,1),coor(i,2),coor(i,3))=1;
        else
            stack(coor(i,1),coor(i,2),coor(i,3))=1;
        end
    end
    clear('coor','i');
end

function stack = fill_stack_2d(coor,nx,ny,valmode)
% fill in a stack [nx ny] with zeros except at positions indicated by coor = [x y I] 
%or [x y z].
%if input positions falls outside of the size of the array, I change the coordinate
%to zero or nx/ny/nz
    if size(coor,2) == 2
        stack = zeros(nx,ny,'uint8');
    elseif size(coor,2) == 3
        stack = zeros(nx,ny,class(coor));
    end

    if ndims(coor)<2
        return
    end  

    coor(:,1) = max(1,min(uint16(ceil(coor(:,1))),nx));
    coor(:,2) = max(1,min(uint16(ceil(coor(:,2))),ny));
    
    for i = 1:size(coor,1)
        if size(coor,2)>=3 && strcmp(valmode,'value')            
            stack(coor(i,1),coor(i,2))=coor(i,3);
        elseif strcmp(valmode,'ones')  
            stack(coor(i,1),coor(i,2))=1;
        else
            stack(coor(i,1),coor(i,2))=1;
        end
    end
    clear('i','coor');
end
