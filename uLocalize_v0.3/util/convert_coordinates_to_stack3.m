function stack = convert_coordinates_to_stack3(varargin)

%convert_coordinates_to_stack([x y z I],nx,ny,nz)
%convert_coordinates_to_stack([x y z],nx,ny,nz)
%convert_coordinates_to_stack([x y z I],nx,ny,nz,'ValMode','ones')
%convert_coordinates_to_stack([x y z I],nx,ny,nz,'ValMode','ones','Size',size)

%convert_coordinates_to_stack([x y I],nx,ny)
%convert_coordinates_to_stack([x y],nx,ny)
%convert_coordinates_to_stack([x y I],nx,ny,'ValMode','ones')
%convert_coordinates_to_stack([x y I],nx,ny,'ValMode','ones','Size',size)

%valmode optional argument
%valmode = 'value' (default): I build a stack with the values corresponding to the last
%column of the array [x y z I] or [x y I]

%valmode = 'ones': I put ones in the stack at the positions indicated by the
%input array [x y z I] or [x y I]

%size : size of the square around the center of the spot 
%I take the closest odd integer from that number

%%
if nargin == 0
    stack = 0;
    return;
end

%% parsing arguments for the dimension of the stack
if size(varargin{1},2) >= 4
    if nargin >= 4
        numdim =3;
        coor = varargin{1};
        nx = varargin{2};
        ny = varargin{3};
        nz = varargin{4};
    else
        stack = 0;
        return;
    end
elseif size(varargin{1},2) == 2
    if nargin >=3
        numdim =2;
        coor = varargin{1};
        nx = varargin{2};
        ny = varargin{3};
    else
        stack = 0;
        return
    end
elseif size(varargin{1},2) == 3
    if nargin <=2
        stack = 0;
        return
    elseif nargin == 3
        numdim = 2;
        coor = varargin{1};
        nx = varargin{2};
        ny = varargin{3};
    elseif nargin >=4
        if isa(varargin{4},'char')
            numdim =2;
            coor = varargin{1};
            nx = varargin{2};
            ny = varargin{3}; 
        else
            numdim = 3;
            coor = varargin{1};
            nx = varargin{2};
            ny = varargin{3};
            nz = varargin{4};
        end
    end
end

%% parsing argument for options
i = 1;
Size = 1;
ValMode = 'ones';
while i<= nargin
    if isa(varargin{i},'char')
        
        if strcmp(varargin{i},'ValMode')
            if i < nargin
                ValMode = varargin{i+1};
                i=i+2;
            end
            
        elseif strcmp(varargin{i},'Size')
            if i < nargin
                Size = varargin{i+1};
                i=i+2;
            end
        else
            i = i+1;
        end
    else 
        i = i+1;
    end
end

if ~isa(Size,'numeric')
    Size = 1;
else
    %the closest odd number (so that the square is centered on the position
    %of the spot
    Size = 2*round(Size/2)+1;
    
end

if ~isa(ValMode,'char')
    ValMode = 'ones';
end

%% generating the stack / image
if numdim == 3
    stack = fill_stack_3d(coor,nx,ny,nz,ValMode,Size);
elseif numdim == 2
     stack = fill_stack_2d(coor,nx,ny,ValMode,Size);
end

end


function stack = fill_stack_3d(coor,nx,ny,nz,ValMode,Size)
% fill in a 3D stack (1:nx,1:ny,1:nz) with zeros except at positions indicated by coor = [x y z I] 
%or [x y z].
%if input positions falls outside of the size of the array, I change the coordinate
%to zero or nx/ny/nz
    
    if size(coor,2) == 3 || strcmp('ValMode','ones')
        stack = zeros(nx,ny,nz,'uint8');
    else
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
        if size(coor,2)>=4 && strcmp(ValMode,'value')            
            stack(:,:,coor(i,3)) = add_spot_to_img(stack(:,:,coor(i,3)),coor(i,1),coor(i,2),coor(i,4),Size);
        elseif strcmp(ValMode,'ones')  
            stack(:,:,coor(i,3)) = add_spot_to_img(stack(:,:,coor(i,3)),coor(i,1),coor(i,2),1,Size);
        else
            stack(:,:,coor(i,3)) = add_spot_to_img(stack(:,:,coor(i,3)),coor(i,1),coor(i,2),1,Size);
        end
    end
    
    if strcmp(ValMode,'ones')  
        if max(stack(:))~=0
            stack = 65535*stack/max(stack(:));
        end
    end
    clear('coor','i');
end

function stack = fill_stack_2d(coor,nx,ny,ValMode,Size)
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
        if size(coor,2)>=3 && strcmp(ValMode,'value')            
            stack = add_spot_to_img(stack,coor(i,1),coor(i,2),coor(i,3),Size);
        elseif strcmp(ValMode,'ones')  
            stack = add_spot_to_img(stack,coor(i,1),coor(i,2),1,Size);
        else
            stack = add_spot_to_img(stack,coor(i,1),coor(i,2),1,Size);
        end
    end
    
    if strcmp(ValMode,'ones')  
        if max(stack(:))~=0
            stack = 65535*stack/max(stack(:));
        end
    end
    clear('i','coor');
end

function img = add_spot_to_img(img,xc,yc,value,Size)

if Size <= 1    
    img(xc,yc) = img(xc,yc)+ value;
else
    [nx,ny] = size(img);
    dx = floor(Size/2);
    x1 = max(1,xc - dx);
    x2 = min(nx,xc + dx);
    y1 = max(1,yc - dx);
    y2 = min(ny,yc + dx);

    img(x1:x2,y1:y2) = img(x1:x2,y1:y2)+ value;
end

end
