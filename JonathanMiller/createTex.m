function createTex(parFile, subject, session)

[pathstr, name, ext] = fileparts(parFile);
[times, type, mode, x, y, angle, speed] = textread(parFile, '%f%s%s%f%f%f%f');
[xStore, yStore] = textread(strcat(parFile,'.str'), '%*s%f%f','headerlines',1);

fid = fopen(strcat(pathstr,'/report.tex'),'w');

fprintf(fid, '\\documentclass{article}\n');
fprintf(fid, '\\usepackage{graphicx}\n');
fprintf(fid, '\\usepackage{color}\n');
fprintf(fid, '\\usepackage[hmargin=3.5cm,vmargin=3cm]{geometry}\n');
fprintf(fid, '\\pagecolor{white}\n');

fprintf(fid, '\\title{BZFlag 2.1 Report}\n');
fprintf(fid, '\\begin{document}\n');

fprintf(fid, '\\maketitle\n');

fprintf(fid, '\\section*{Session Information}\n');

fprintf(fid, 'Subject: %s\\\\\n',subject);
fprintf(fid, 'Session: %s\\\\\n',session(end));

fprintf(fid, '\\section*{Basic Statistics}\n');

[trainingTrials, meanTraining, seekTrials, meanSeek] = basicStats(times, type, mode);

fprintf(fid, '\\subsection*{Training Phase}\n');
fprintf(fid, 'Trials completed: %d\\\\\n', trainingTrials);
fprintf(fid, 'Mean trial time: %0.2fs\\\\\n', meanTraining);

fprintf(fid, '\\subsection*{Seek Phase}\n');
fprintf(fid, 'Trials completed: %d\\\\\n', seekTrials);
fprintf(fid, 'Mean trial time: %0.2fs\\\\\n', meanSeek);

fprintf(fid, '\\section*{Occupancy Plot}\n');

occFile = occupancy(times, type, mode, x, y, speed, pathstr);

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\includegraphics[width=100mm]{');
fprintf(fid, occFile);
fprintf(fid, '}\\\\\n');
fprintf(fid, '\\end{center}\n');

fprintf(fid, '\\section*{Learning Curve by Latency per Trial}\n');

curve = latencyLearningCurve(times, type, mode);

[curveFile, histFile] = combineCurves(curve, 'latencyLearningCurve', 'Trial', 'Latency (s)', pathstr);

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\includegraphics[width=100mm]{');
fprintf(fid, curveFile);
fprintf(fid, '}\\\\\n');
fprintf(fid, '\\end{center}\n');


fprintf(fid, '\\section*{Learning Curve by Latency per Round}\n');

curve = blockedLatencyLearningCurve(times, type, mode);
[curveFile,histFile] = combineCurves(curve, 'blockedLatencyLearningCurve', 'Round', 'Latency (s)', pathstr);

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\includegraphics[width=100mm]{');
fprintf(fid, curveFile);
fprintf(fid, '}\\\\\n');
fprintf(fid, '\\end{center}\n');



fprintf(fid, '\\section*{Learning Curve by Excess Distance per Trial}\n');

curve = errorLearningCurve(times, type, mode, x, y);
[curveFile,histFile] = combineCurves(curve, 'errorLearningCurve', 'Trial', 'Excess Distance (VR Units)', pathstr);

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\includegraphics[width=100mm]{');
fprintf(fid, curveFile);
fprintf(fid, '}\\\\\n');
fprintf(fid, '\\end{center}\n');

fprintf(fid, '\\subsection*{Histogram of Excess Distance}\n');

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\includegraphics[width=100mm]{');
fprintf(fid, histFile);
fprintf(fid, '}\\\\\n');
fprintf(fid, '\\end{center}\n');



fprintf(fid, '\\section*{Learning Curve by Excess Distance per Round}\n');

curve = blockedErrorLearningCurve(times, type, mode, x, y);
[curveFile, histFile] = combineCurves(curve, 'blockedErrorLearningCurve', 'Round', 'Excess Distance (VR Units)', pathstr);
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\includegraphics[width=100mm]{');
fprintf(fid, curveFile);
fprintf(fid, '}\\\\\n');
fprintf(fid, '\\end{center}\n');


fprintf(fid, '\\end{document}\n');

fclose(fid);

quit;
