%% it collects the data, for shape and appearance and produce 3 scales and
%% 4 rotations for shape and apperance.
function [appMatrix, shMatrix]= fbSrPerturbScaleAndRotation( imageTensor, complexCoords, showIt )


[ numLm numSamp ] = size( complexCoords );



%shape part
% scales=[.8 1 1.2];
angles=[ -8 0 8 ];

scales=[ 1 ];
% transX=[0];
% transY=[0];

% initialize
shMatrix=zeros(numLm,numSamp*length(scales)*length(angles));
appMatrix = zeros( size(imageTensor,1), size(imageTensor,2), ...
                   numSamp*length(scales)*length(angles));
% trMatrix=zeros(2,numSamp);
% sc_angMatrix=zeros(2,numSamp);


%% appearance part
% rszWind=150;
% appMatrix=zeros([rszWind^2,numSamp*3*5],'single');


idx=1;
for k1=1:numSamp,
    
    %im1=imread([imFolder maList{k1}(1:end-4) '.bmp']);    
    %load([maFolder maList{k1}]);
    im1scl = imageTensor(:,:,k1);
    faceCoordinatesScl = [ real( complexCoords(:,k1)) imag( complexCoords(:,k1)) ];
    
    
    for k2=scales,

%         cropConstant=k2;
%         maxX=max(faceCoordinates(:,1))+cropConstant;
%         minX=min(faceCoordinates(:,1))-cropConstant;
%         maxY=max(faceCoordinates(:,2))+cropConstant;
%         minY=min(faceCoordinates(:,2))-cropConstant;
%         
%         im1scl=imcrop(im1,[minX minY maxX-minX maxY-minY]);
%         faceCoordinatesScl(:,1)=faceCoordinates(:,1)-minX;
%         faceCoordinatesScl(:,2)=faceCoordinates(:,2)-minY;
% 
%         [height width l]=size(im1scl);
        
        %-------------------------
        % added by Sam
        
        %show before scaling
        % imagesc( im1scl), colormap gray, hold on;
        % plot( faceCoordinatesScl(:,1), faceCoordinatesScl(:,2), 'g*');
        % pause;
        
        [height width ]=size(im1scl); 
        
        % rescale
        im1scl = imresize( im1scl, k2);   
        afterScaleSize = size( im1scl);
        faceCoordinatesScl = k2.*faceCoordinatesScl;
        % crop it, keep original input image size
        [ im1scl complexCoordsCropped ] = cropFiducial( im1scl, [faceCoordinatesScl(:,1) + 1i*faceCoordinatesScl(:,2)], ...
                                       1:size(faceCoordinatesScl,1), [height width ], 0, ...
                                       [], 0, .5*[afterScaleSize(1) afterScaleSize(2)] );
        faceCoordinatesScl = [real(complexCoordsCropped) imag(complexCoordsCropped) ];
                                                 
        % show after scaling
        % imagesc( im1scl), colormap gray, hold on;
        % plot( faceCoordinatesScl(:,1), faceCoordinatesScl(:,2), 'g*');
        % pause;
        %---------------------------------------
        
        for k3=angles,
%             for k4=1:length(transX),
%                 for k5=1:length(transY),
                    im1Temp=im1scl;
                    lm=faceCoordinatesScl;
                    im1Temp=imrotate(im1Temp,k3,'bilinear','crop');                    
                    lm(:,1)=lm(:,1)-width/2;
                    lm(:,2)=lm(:,2)-height/2;
                    lm=lm*...
                        [cos(-k3*pi/180) sin(-k3*pi/180);-sin(-k3*pi/180) cos(-k3*pi/180)];                    
                    lm(:,1)=lm(:,1)+width/2;
                    lm(:,2)=lm(:,2)+height/2;  
                    
                    % show after rotation
                    if showIt
                        clf(gcf)
                        imagesc( im1Temp), colormap gray, hold on;
                        plot( lm(:,1), lm(:,2), 'g*');
                        pause(.02);
                    end        
                    
%                     im1Temp=padarray(im1Temp,[transY(k5) transX(k4)],'pre');
%                     [heightTmp widthTmp l]=size(im1Temp);
%                     hdel=round(abs(height-heightTmp)/2)+1;
%                     wdel=round(abs(width-widthTmp)/2)+1;
%                     
%                     lm(:,1)=lm(:,1)+transX(k4)-wdel+1;
%                     lm(:,2)=lm(:,2)+transY(k5)-hdel+1;
%                     
%                     
%                     im1Temp=im1Temp(hdel:end-hdel,wdel:end-wdel,:);
%                     [h2 w2 l]=size(im1Temp);
%                     
%                     
%                     lm(:,1)=lm(:,1)*rszWind/w2;
%                     lm(:,2)=lm(:,2)*rszWind/h2;
%                     im1Temp=rgb2gray(imresize(im1Temp,[rszWind rszWind]));
% 
%                     imagesc(im1Temp);hold on;colormap gray
%                     plot(lm(:,1),lm(:,2),'.')
%                     hold off;
%                     im1Data=reshape(im1Temp,rszWind^2,1);

                    appMatrix(:,:,idx)=im1Temp; % single(single(im1Data)./norm(single(im1Data)));
                    shMatrix(:,idx)=lm(:,1) + 1i*lm(:,2);
                    
                    idx=idx+1;
                    
%                 end
%             end
%           imagesc(im1Temp(:,:,1));
        end
    end
end



