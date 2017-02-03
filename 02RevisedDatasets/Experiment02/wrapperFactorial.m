%resolutionPop  = 5;
%resolutionRho  = 30;
repetitions = 20;
stepSize=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=linspace(.98,1,stepSize)
y=logspace(-2,-7,stepSize)
rhoVector=x.*y
eVector=x.*(1-y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sVector     = 0;                                            %linspace(0,1,resolution);
%rhoVector   = [0 logspace(-7,-3,resolutionRho)];            %linspace(0,1,resolution);
%eVector     = .98;                                          %linspace(0,1,resolution);
MeqValue    = 10000%[1000 5000 10000 25000 50000 100000]%linspace(1000,100000,resolutionPop);           %1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the number of iterations to run
iter=0;
    for j=1:length(eVector)
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
    for j=1:length(eVector)
        for i=1:length(rhoVector)
            parfor repetition = 1:repetitions
                iterateOneSimStep(rhoVector(i),eVector(j),sVector,MeqValue,int2str(randi([0 1000000000],1)),'./Data2/');
            end
            iter=iter+repetitions;
            iter/totalIterations
        end
    end
%end  
%%

