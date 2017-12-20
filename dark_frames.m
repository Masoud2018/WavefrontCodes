

function [time1,time2] = dark_frames(Undulator,Regime)


if (Undulator == 0 && Regime == 0 ) %A1 mode
    
  % start time
  time1 = '2015-06-28 23:29:36';
  % end time 
  time2 = '2015-06-28 23:48:16';
    
else
   if (Undulator == 1 && Regime == 0 ) %A2 mode
        
        % start time
        time1 = '2015-06-28 18:29:22';
        % end time 
        time2 = '2015-06-28 18:35:49';
   else  
      if (Undulator == 0 && Regime == 1 ) % A3 mode
            
            % start time
            time1 = '2015-06-29 02:54:04';
            % end time 
            time2 = '2015-06-29 03:02:42';
      else  
          if (Undulator == 1 && Regime == 1 ) % A4 mode
            
              % start time
              time1 = '2015-06-22 17:12:51';
              % end time 
              time2 = '2015-06-22 17:16:53';
              
          else
              fprintf('numbers are not correct')
              
          end      
      end
   end
end
end