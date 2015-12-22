load SaccDistribution


qads = {[1:10], [11:22], 23:32}
ctr = 0;
for i = 1:3
     for j = 1:3
          ctr = ctr + 1;
          redSaccDist{i,j} = zeros(32, 32);
          for x=qads{i}
               for y = qads{j}
                    redSaccDist{i,j} = redSaccDist{i,j}+reshape(SaccByPos(x,y,:,:), [32 32]);
               end
          end
          subplot(3,3,ctr)
          imshow(redSaccDist{i,j}, []);
     end
end
colormap('hot')

subplot(3,3,1);
title('saccades from top left region')
subplot(3,3,2);
title('saccades from top region')
subplot(3,3,3);
title('saccades from top right region')
subplot(3,3,4);
title('saccades from left region')
subplot(3,3,5);
title('saccades from centre region')
subplot(3,3,6);
title('saccades from right region')
subplot(3,3,7);
title('saccades from lower left region')
subplot(3,3,8);
title('saccades from lower region')
subplot(3,3,9);
title('saccades from lower right region')
set(gcf, 'Color', 'w');
export_fig saccDistExample.pdf