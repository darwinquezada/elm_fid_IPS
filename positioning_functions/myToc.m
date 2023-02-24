function [elapsedTime] = myToc(lastTicMeasuredTime)
%Joaquín Torres-Sospedra
  elapsedTime = cputime - lastTicMeasuredTime;
  %fprintf('Elapsed cpu time is %f seconds.\n',elapsedTime);
end