%% Restart function

function []=defaultconfiguration(app)
    % Clear XTick and YTick for Ax2 and Ax3
    app.Ax2.XTick=[];
    app.Ax2.YTick=[];
    app.Ax3.XTick=[];
    app.Ax3.YTick=[];
    
    %Enable panels
    panelsToEnable = [app.LoadImagePanel, app.LoadImagePanel_2];
    set(findall(panelsToEnable, '-property', 'enable'), 'enable', 'on');

    
    % Load image and source files
    files = dir(fullfile(app.imfolder, '*.tif'));
    ListItems = {files.name};

    files2 = dir(fullfile(app.mainfolder, '*.xlsx'));
    Sourcefiles = {files2.name};
    
    % Update dropdown items
    app.ImageDropDown.Items=ListItems;
    app.MaskDropDown.Items=ListItems;
    app.ImageDropDown_2.Items=ListItems;
    app.ImageDropDown_3.Items=ListItems;
    app.ImageDropDown_4.Items=ListItems;
    app.ImageDropDown_5.Items=ListItems;
    app.MaskDropDown_2.Items=ListItems;
    app.SourceDropDown.Items=Sourcefiles;
    
    % Select default tab
    app.TabGroup.SelectedTab = app.OneSimulationTab;
    
end