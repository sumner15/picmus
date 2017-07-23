function out = picmus_path()
%PICMUS_PATH Returns the path of the picmus ToolBox
%   Usage: out = picmus_path()
    
    [out,name,ext] = fileparts(which('picmus_path'));
end
