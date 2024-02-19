function writeGIF(dataset)
    filename = 'movie.gif';

    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    tempRes = 0.035;
    frames = size(dataset,3);
    lowerLim = 0; %window level
    upperLim = 1800;
    for i=1:frames
        imshow(dataset(:,:,i),[lowerLim upperLim]); drawnow;
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        if i == 1 
          imwrite(imind,cm,filename,'gif','DelayTime',tempRes,'Loopcount',inf); 
        else 
          imwrite(imind,cm,filename,'gif','DelayTime',tempRes,'WriteMode','append'); 
        end 
    end
end 

