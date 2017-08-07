%
% plotbumptraj.m
%

clear stdpatch
vertx=[[0:dt:T-200]'; [T-200:-dt:0]']; vertx(:,2)=[sqrt(cmv)'; -sqrt(cmv(end:-1:1))';];
stdpatch.Vertices=vertx;
stdpatch.Faces=[1:length(vertx)];
stdpatch.FaceColor=[250 128 114]/255;
stdpatch.EdgeColor=[250 128 114]/255;
figure(2), hold on, patch(stdpatch);

for j=1:16,
    hold on, plot([0:dt:T-200],dcen(j,:),'b','linewidth',1);
end

set(gca,'xtick',[0:100:600]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);