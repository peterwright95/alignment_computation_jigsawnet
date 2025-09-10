function save_fig( filename )
    old_fig_size = get(gcf,'Position');
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf, 'Color', 'white');
    out_ext='jpg';
    [outdir,output_name,~]=fileparts(ufilename(filename));
    export_fig( gcf,fullfile(outdir,output_name),['-' out_ext],'-nocrop','-q75'); 
    set(gcf, 'Color', 'default');
    set(gcf, 'Position', old_fig_size);
end

