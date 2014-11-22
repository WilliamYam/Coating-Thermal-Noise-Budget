classdef MatlabPlotterFactory < handle
    
    properties
    end
    
    methods
        function plotter = getPlotter(self, noiseModel)
            plotter = NoisePlotter(noiseModel);
            plotter.epilog{end+1} = @self.stripLinks;
        end
        
        function stripLinks(~, noisePlotter, ~)
            lineObjs = noisePlotter.handles.ln;
            for n = 1:numel(lineObjs)
                str = get(lineObjs(n), 'DisplayName');
                % strip href and hyperlink links (\href, \hyperlink)
                % assume no latex inside the link text
                str = regexprep(str, '\\href{[^}]*}{([^\\]*)}', '$1');
                str = regexprep(str, '\\hyperlink{[^}]*}{([^\\]*)}', '$1');
                set(lineObjs(n), 'DisplayName', str);
            end
            textObjs = findall(noisePlotter.handles.fg, 'Type', 'text');
            for n = 1:numel(textObjs)
                str = get(textObjs(n), 'String');
                % string with line breaks is returned as a 2D array
                % convert such strings to 1D array for regexprep
                strL = cellstr(str);
                str = [sprintf('%s\n', strL{1:end-1}), strL{end}];
                % strip hyperlink target (\hypertarget)
                % assume target wraps the entire string, and may contain
                % links inside the target text
                str = regexprep(str, '\\hypertarget{[^}]*}{(.*)}', '$1');
                % strip href and hyperlink links (\href, \hyperlink)
                % assume no latex inside the link text
                str = regexprep(str, '\\href{[^}]*}{([^\\]*)}', '$1');
                str = regexprep(str, '\\hyperlink{[^}]*}{([^\\]*)}', '$1');
                set(textObjs(n), 'String', str);
            end
        end
        
        function finalize(~)
        end
        
        function render(~)
        end
        
        function cleanup(~)
        end
    end
    
end

