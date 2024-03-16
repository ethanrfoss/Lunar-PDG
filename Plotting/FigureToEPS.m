function FigureToEPS(SaveName)

currentFile = mfilename('fullpath');
currentDirectory = fileparts(currentFile);

print([currentDirectory '/SavedImages/' SaveName],'-depsc2');

end