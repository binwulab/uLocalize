function img=mytiffread(filename, index)
if exist(filename,'file') ~= 2
    disp(['mytifread: file does not exist: ', filename]);
    img=[];
    return;
end
imf=imfinfo(filename);
n=numel(imf);
if nargin<2
    index=1:n;
end
if find(ismember(index, 1:n)==0)
    disp(['mytifread: index exceed the bound: ', filename]);
    img=[];
    return;
end 
for i=1:numel(index)
    if i==1
        img=imread(filename,index(1));
        img=repmat(img,1,1, numel(index));
    else
        img(:,:,i)=imread(filename,index(i));
    end
end
end