resolutionPop  = 5;
resolutionRho  = 30;
repetitions = 128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sVector     = 0;                                            %linspace(0,1,resolution);
rhoVector   = .01;                        %linspace(0,1,resolution);
eVector     = .98;                                          %linspace(0,1,resolution);
MeqValue    = 10000;                         %linspace(1000,100000,resolutionPop);           %1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the number of iterations to run
iter=0;
    for j=1:length(MeqValue)
        for i=1:length(rhoVector)
            for repetition = 1:repetitions
               iter=iter+1;
            end
        end
    end
totalIterations=iter
%% 
%Run iterations in parallel
iter=0;
%for k=1:resolution
    %for j=1:length(MeqValue)
        %for i=1:length(rhoVector)
            parfor repetition = 1:repetitions
                iterateOneSimStep(rhoVector,eVector,sVector,MeqValue,int2str(randi([0 1000000000],1)),'./DataB/');
            end
            %iter=iter+repetitions;
            %iter/totalIterations
        %end
    %end
%end  