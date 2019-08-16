function [coeffDic, TDic, RDic] = loadParameters( xmlFileName )
% Load parameters of different types of SMA from xml file
% -- input
% xmlFileName -- file name of xml file
% -- output
% coeffDic -- coefficient parameters dictionary
% TDic     -- transformation temperatures parameter dictionary

if (nargin < 1)
    xmlFileName = 'default.xml';
end

try
    evalin('base', 'loadParamsFlag;');
catch
    assignin('base', 'loadParamsFlag', 1);
    settings = xmlread(xmlFileName, 'r');
    parametersNode = settings.getDocumentElement();

    % cosntant mechanical parameters
    % coefficients
    D     = str2double(parametersNode.getElementsByTagName('YoungsModulusInMPa').item(0).getTextContent());
    THETA = str2double(parametersNode.getElementsByTagName('ThermoelasticTensorInMPaPerDegreeCelsius').item(0).getTextContent());
    OMEGA = str2double(parametersNode.getElementsByTagName('TransformationTensorInMPa').item(0).getTextContent());

    % transformation temperatures
    M_f = str2double(parametersNode.getElementsByTagName('MartensiteFinishInDegreeCelsius').item(0).getTextContent());
    M_s = str2double(parametersNode.getElementsByTagName('MartensiteStartInDegreeCelsius').item(0).getTextContent());
    A_s = str2double(parametersNode.getElementsByTagName('AusteniteStartInDegreeCelsius').item(0).getTextContent());
    A_f = str2double(parametersNode.getElementsByTagName('AusteniteFinishInDegreeCelsius').item(0).getTextContent());

    % transition temperature stress relationship
    C_A = str2double(parametersNode.getElementsByTagName('C_A').item(0).getTextContent());
    C_M = str2double(parametersNode.getElementsByTagName('C_M').item(0).getTextContent());

    coeffDic = struct('YoungsModulus', D, 'ThermoelasticTensor', THETA, 'TransformationTensor', OMEGA);
    TDic     = struct('MartensiteFinish', M_f, 'MartensiteStart', M_s, 'AusteniteStart', A_s, 'AusteniteFinish', A_f);
    RDic     = struct('C_A', C_A, 'C_M', C_M);

    assignin('base', 'coeffDic', coeffDic);
    assignin('base', 'TDic', TDic);
    assignin('base', 'RDic', RDic);
    
    return;
end

coeffDic = evalin('base', 'coeffDic');
TDic     = evalin('base', 'TDic');
RDic     = evalin('base', 'RDic');

end