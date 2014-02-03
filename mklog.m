function mklog(messageToLog)
display(sprintf('%s - %s',datestr(now),messageToLog));
f=fopen('log.txt','a');
fprintf(f,sprintf('%s - %s\r\n',datestr(now),messageToLog));
fclose(f);