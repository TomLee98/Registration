function status = ocreg(movsrc_, movdst_, movtmpl_, regfrs_, regopt_)
%OCREG This function is one channels registration algorithm, which
%corresponds to the caller RegisterWorker
% Input:
%   - movsrc_: 1-by-1 regmov object, source registration movie
%   - movdst_: 1-by-1 regmov object, destination registration movie
%   - movtmpl_:1-by-1 uint16 array, the registration template volume
%   - regfrs_: 1-by-n positive integer numeric array, indicating the frames
%              need to be aligned
%   - regopt_: 1-by-1 regopt object, with registration options
% Output:
%   - status: 1-by-1 double, 0 for normal exit, other code for bad exit
%
%   see also: regmov, regtmpl, imregtform, imregcorr, imwarp, imregdemons, 
%   imhistmatchn, imref3d

arguments
    movsrc_ (1,1)   regmov
    movdst_ (1,1)   regmov
    movtmpl_(:,:,:) uint16
    regfrs_ (1,:)   double {mustBePositive, mustBeInteger}
    regopt_ (1,1)   struct {mustBeRegistrationOption}
end
end

function mustBeRegistrationOption(A)
VALID_FIELD = ["Mode", "RegModal", "TformType", "MaxStep", "MinStep", "MaxIterN", ...
    "IterCoeff", "VPL", "Interp", "CoRegC", "DS", "SC", "FC", "Hardware"];

if ~all(ismember(fieldnames(A), VALID_FIELD))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

% Escape value checking for tcreg calling faster

end