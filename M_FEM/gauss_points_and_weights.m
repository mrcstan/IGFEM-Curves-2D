%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 11/3/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns guss points and weights for triangular and
% quadrilateral elements
% isTriangular: true - the integration domain is triangular or tetrahedral 
% npt: total number of gauss points for triangular or tetrahedral
%      integration. number of gauss points in one direction for quadrilateral 
%      or hexahedral integration
% format: 'combined': combine the points and the weights
%         'separate': seperate the points and the weights
% OUTPUT:
%   if format == 'combined'
%       gauss = [x11,x12,...,x1n;x21,x22,...,x2n;
%               ...;xm1,xm2,...,xmn;w1,w2,...,wn];
%               where m is the number of dimensions and n is the number of
%               gauss points
%   if format == 'separate'
%       gauss.pt = [x11,x12,...,x1n;x21,x22,...,x2n;
%               ...;xm1,xm2,...,xmn];
%       gauss.weight = [w1,w2,...,wn];   
function gauss = gauss_points_and_weights(isTriangular,npt,dimension,format) 
if (nargin < 4)
    formatID = 0;
else
    if (strcmpi(format,'separate'))
        formatID = 0;
    elseif (strcmpi(format,'combined'))
        formatID = 1;
    else
        error('format not recognized')
    end
end
if (isTriangular)
    switch dimension
        case 1
             [pt,weight]=lgwt(npt,0,1);
             pt = pt';
             weight = weight';
        case 2
            switch npt
                case 3 % error is of order O(h^3)
                    pt = [1./6.,1./6.;
                                   2./3.,1./6.;
                                   1./6.,2./3.]';
                    weight = [1./6., 1./6., 1./6.];
                case 7 % error is of order O(h^6)
                    r1=0.1012865073235;
                    r2=0.7974269853531;
                    r4=0.4701420641051;
                    r6=0.0597158717898;
                    r7=1./3.;
                    pt = [r1,r1;
                         r2,r1;
                         r1,r2;
                         r4,r6;
                         r4,r4;
                         r6,r4;
                         r7,r7]';
                    w1=0.1259391805448;
                    w4=0.1323941527885;
                    w7=0.225;
                    weight = 0.5*[w1,w1,w1,w4,w4,w4,w7];
                case 13 % error is of order O(h^8)
                    ra1 = 0.333333333333333;
                    rb1 = 0.479308067841920;
                    rb2 = 0.260345966079040;
                    rc1 = 0.869739794195568;
                    rc2 = 0.065130102902216;
                    rd1 = 0.048690315425316;
                    rd2 = 0.312865496004874;
                    rd3 = 0.638444188569810;
                    pt = [ra1,ra1;
                          rb1,rb2;
                          rb2,rb1;
                          rb2,rb2;
                          rc1,rc2;
                          rc2,rc1;
                          rc2,rc2;
                          rd1,rd2;
                          rd1,rd3;
                          rd2,rd3;
                          rd2,rd1;
                          rd3,rd1;
                          rd3,rd2]';
                    wa = -0.149570044467682;
                    wb = 0.175615257433208;
                    wc = 0.053347235608838;
                    wd = 0.077113760890257;
                    weight = 0.5*[wa,wb,wb,wb,wc,wc,wc,wd,wd,wd,wd,wd,wd];
				case 16 % error is of order O(h^9)
					pt = [3.3333333333333300E-01	3.3333333333333300E-01
                        8.1414823414554000E-02	4.5929258829272300E-01
                        4.5929258829272300E-01	4.5929258829272300E-01
                        4.5929258829272300E-01	8.1414823414554000E-02
                        6.5886138449648000E-01	1.7056930775175900E-01
                        1.7056930775175900E-01	1.7056930775175900E-01
                        1.7056930775175900E-01	6.5886138449648000E-01
                        8.9890554336593700E-01	5.0547228317030900E-02
                        5.0547228317030900E-02	5.0547228317030900E-02
                        5.0547228317030900E-02	8.9890554336593700E-01
                        8.3947774099580000E-03	2.6311282963463800E-01
                        2.6311282963463800E-01	7.2849239295540400E-01
                        7.2849239295540400E-01	8.3947774099580000E-03
                        2.6311282963463800E-01	8.3947774099580000E-03
                        7.2849239295540400E-01	2.6311282963463800E-01
                        8.3947774099580000E-03	7.2849239295540400E-01]';
					weight = 0.5*[1.4431560767778700E-01
                                9.5091634267284900E-02
                                9.5091634267284900E-02
                                9.5091634267284900E-02
                                1.0321737053471700E-01
                                1.0321737053471700E-01
                                1.0321737053471700E-01
                                3.2458497623198000E-02
                                3.2458497623198000E-02
                                3.2458497623198000E-02
                                2.7230314174435000E-02
                                2.7230314174435000E-02
                                2.7230314174435000E-02
                                2.7230314174435000E-02
                                2.7230314174435000E-02
                                2.7230314174435000E-02]';
					
                case 25 % error is of order O(h^11)
                    ra1 = 0.333333333333333;
                    rb1 = 0.028844733232685;
                    rb2 = 0.485577633383657;
                    rc1 = 0.781036849029926;
                    rc2 = 0.109481575485037;
                    rd1 = 0.141707219414880;
                    rd2 = 0.307939838764121;
                    rd3 = 0.550352941820999;
                    re1 = 0.025003534762686;
                    re2 = 0.246672560639903;
                    re3 = 0.728323904597411;
                    rf1 = 0.009540815400299;
                    rf2 = 0.066803251012200;
                    rf3 = 0.923655933587500;
                    wa = 0.090817990382754;
                    wb = 0.036725957756467;
                    wc = 0.045321059435528;
                    wd = 0.072757916845420;
                    we = 0.028327242531057;
                    wf = 0.009421666963733;
                    pt = [ra1,ra1;
                          rb1,rb2;
                          rb2,rb1;
                          rb2,rb2;
                          rc1,rc2;
                          rc2,rc1;
                          rc2,rc2;
                          rd1,rd2;
                          rd1,rd3;
                          rd2,rd3;
                          rd2,rd1;
                          rd3,rd1;
                          rd3,rd2;
                          re1,re2;
                          re1,re3;
                          re2,re3;
                          re2,re1;
                          re3,re1;
                          re3,re2;
                          rf1,rf2;
                          rf1,rf3;
                          rf2,rf3;
                          rf2,rf1;
                          rf3,rf1;
                          rf3,rf2]';
                    weight = 0.5*[wa,wb,wb,wb,wc,wc,wc,wd,wd,wd,wd,wd,wd,...
                                  we,we,we,we,we,we,wf,wf,wf,wf,wf,wf];
				case 33 % error is of order O(h^13)
					pt = [2.3565220452390000E-02	4.8821738977380500E-01
							4.8821738977380500E-01	4.8821738977380500E-01
							4.8821738977380500E-01	2.3565220452389900E-02
							1.2055121541107800E-01	4.3972439229445900E-01
							4.3972439229445900E-01	4.3972439229445900E-01
							4.3972439229445900E-01	1.2055121541107800E-01
							4.5757922997576800E-01	2.7121038501211500E-01
							2.7121038501211500E-01	2.7121038501211500E-01
							2.7121038501211500E-01	4.5757922997576800E-01
							7.4484770891682800E-01	1.2757614554158600E-01
							1.2757614554158600E-01	1.2757614554158600E-01
							1.2757614554158600E-01	7.4484770891682800E-01
							9.5736529909357900E-01	2.1317350453210000E-02
							2.1317350453210000E-02	2.1317350453210000E-02
							2.1317350453210000E-02	9.5736529909357900E-01
							1.1534349453469800E-01	2.7571326968551400E-01
							2.7571326968551400E-01	6.0894323577978800E-01
							6.0894323577978800E-01	1.1534349453469800E-01
							2.7571326968551400E-01	1.1534349453469800E-01
							6.0894323577978800E-01	2.7571326968551400E-01
							1.1534349453469800E-01	6.0894323577978800E-01
							2.2838332222257000E-02	2.8132558098993900E-01
							2.8132558098993900E-01	6.9583608678780200E-01
							6.9583608678780200E-01	2.2838332222257000E-02
							2.8132558098993900E-01	2.2838332222257000E-02
							6.9583608678780200E-01	2.8132558098993900E-01
							2.2838332222257000E-02	6.9583608678780200E-01
							2.5734050548330000E-02	1.1625191590759600E-01
							1.1625191590759600E-01	8.5801403354407300E-01
							8.5801403354407300E-01	2.5734050548330000E-02
							1.1625191590759600E-01	2.5734050548330000E-02
							8.5801403354407300E-01	1.1625191590759600E-01
							2.5734050548330000E-02	8.5801403354407300E-01]';
					weight = 0.5*[2.5731066440455000E-02
								2.5731066440454900E-02
								2.5731066440454900E-02
								4.3692544538038000E-02
								4.3692544538038000E-02
								4.3692544538038000E-02
								6.2858224217885000E-02
								6.2858224217885000E-02
								6.2858224217885000E-02
								3.4796112930709000E-02
								3.4796112930709000E-02
								3.4796112930709000E-02
								6.1662610515590000E-03
								6.1662610515590000E-03
								6.1662610515590000E-03
								4.0371557766381000E-02
								4.0371557766381000E-02
								4.0371557766381000E-02
								4.0371557766381000E-02
								4.0371557766381000E-02
								4.0371557766381000E-02
								2.2356773202303000E-02
								2.2356773202303000E-02
								2.2356773202303000E-02
								2.2356773202303000E-02
								2.2356773202303000E-02
								2.2356773202303000E-02
								1.7316231108659000E-02
								1.7316231108659000E-02
								1.7316231108659000E-02
								1.7316231108659000E-02
								1.7316231108659000E-02
								1.7316231108659000E-02]';
				case 61 % error is of order O(h^18)
					pt =   [3.3333333333333300E-01	3.3333333333333300E-01
							5.6589188864520000E-03	4.9717054055677400E-01
							4.9717054055677400E-01	4.9717054055677400E-01
							4.9717054055677400E-01	5.6589188864520000E-03
							3.5647354750751000E-02	4.8217632262462500E-01
							4.8217632262462500E-01	4.8217632262462500E-01
							4.8217632262462500E-01	3.5647354750751000E-02
							9.9520061958437000E-02	4.5023996902078200E-01
							4.5023996902078200E-01	4.5023996902078200E-01
							4.5023996902078200E-01	9.9520061958437000E-02
							1.9946752124520600E-01	4.0026623937739700E-01
							4.0026623937739700E-01	4.0026623937739700E-01
							4.0026623937739700E-01	1.9946752124520600E-01
							4.9571746405809500E-01	2.5214126797095200E-01
							2.5214126797095200E-01	2.5214126797095200E-01
							2.5214126797095200E-01	4.9571746405809500E-01
							6.7590599068307700E-01	1.6204700465846100E-01
							1.6204700465846100E-01	1.6204700465846100E-01
							1.6204700465846100E-01	6.7590599068307700E-01
							8.4824823547850800E-01	7.5875882260746000E-02
							7.5875882260746000E-02	7.5875882260746000E-02
							7.5875882260746000E-02	8.4824823547850800E-01
							9.6869054606435600E-01	1.5654726967822000E-02
							1.5654726967822000E-02	1.5654726967822000E-02
							1.5654726967822000E-02	9.6869054606435600E-01
							1.0186928826919000E-02	3.3431986736365700E-01
							3.3431986736365700E-01	6.5549320380942300E-01
							6.5549320380942300E-01	1.0186928826919000E-02
							3.3431986736365700E-01	1.0186928826919000E-02
							6.5549320380942300E-01	3.3431986736365700E-01
							1.0186928826919000E-02	6.5549320380942300E-01
							1.3544087167103600E-01	2.9222153779694300E-01
							2.9222153779694300E-01	5.7233759053202000E-01
							5.7233759053202000E-01	1.3544087167103600E-01
							2.9222153779694300E-01	1.3544087167103600E-01
							5.7233759053202000E-01	2.9222153779694300E-01
							1.3544087167103600E-01	5.7233759053202000E-01
							5.4423924290583000E-02	3.1957488542319000E-01
							3.1957488542319000E-01	6.2600119028622800E-01
							6.2600119028622800E-01	5.4423924290583000E-02
							3.1957488542319000E-01	5.4423924290583000E-02
							6.2600119028622800E-01	3.1957488542319000E-01
							5.4423924290583000E-02	6.2600119028622800E-01
							1.2868560833637000E-02	1.9070422419229100E-01
							1.9070422419229100E-01	7.9642721497407000E-01
							7.9642721497407000E-01	1.2868560833637000E-02
							1.9070422419229100E-01	1.2868560833637000E-02
							7.9642721497407000E-01	1.9070422419229100E-01
							1.2868560833637000E-02	7.9642721497407000E-01
							6.7165782413524000E-02	1.8048321164874600E-01
							1.8048321164874600E-01	7.5235100593772800E-01
							7.5235100593772800E-01	6.7165782413524000E-02
							1.8048321164874600E-01	6.7165782413524000E-02
							7.5235100593772800E-01	1.8048321164874600E-01
							6.7165782413524000E-02	7.5235100593772800E-01
							1.4663182224828000E-02	8.0711313679563900E-02
							8.0711313679563900E-02	9.0462550409560800E-01
							9.0462550409560800E-01	1.4663182224828000E-02
							8.0711313679563900E-02	1.4663182224828000E-02
							9.0462550409560800E-01	8.0711313679563900E-02
							1.4663182224828000E-02	9.0462550409560800E-01]';
					weight = 0.5*[3.3437199290803000E-02
								5.0934154405070000E-03
								5.0934154405070000E-03
								5.0934154405070000E-03
								1.4670864527638000E-02
								1.4670864527638000E-02
								1.4670864527638000E-02
								2.4350878353672000E-02
								2.4350878353672000E-02
								2.4350878353672000E-02
								3.1107550868968900E-02
								3.1107550868968900E-02
								3.1107550868968900E-02
								3.1257111218620000E-02
								3.1257111218620000E-02
								3.1257111218620000E-02
								2.4815654339665000E-02
								2.4815654339665000E-02
								2.4815654339665000E-02
								1.4056073070557000E-02
								1.4056073070557000E-02
								1.4056073070557000E-02
								3.1946761737789900E-03
								3.1946761737789900E-03
								3.1946761737789900E-03
								8.1196553189930000E-03
								8.1196553189930000E-03
								8.1196553189930000E-03
								8.1196553189930000E-03
								8.1196553189930000E-03
								8.1196553189930000E-03
								2.6805742283163000E-02
								2.6805742283163000E-02
								2.6805742283163000E-02
								2.6805742283163000E-02
								2.6805742283163000E-02
								2.6805742283163000E-02
								1.8459993210822000E-02
								1.8459993210822000E-02
								1.8459993210822000E-02
								1.8459993210822000E-02
								1.8459993210822000E-02
								1.8459993210822000E-02
								8.4768685343280000E-03
								8.4768685343280000E-03
								8.4768685343280000E-03
								8.4768685343280000E-03
								8.4768685343280000E-03
								8.4768685343280000E-03
								1.8292796770024900E-02
								1.8292796770024900E-02
								1.8292796770024900E-02
								1.8292796770024900E-02
								1.8292796770024900E-02
								1.8292796770024900E-02
								6.6656320041650000E-03
								6.6656320041650000E-03
								6.6656320041650000E-03
								6.6656320041650000E-03
								6.6656320041650000E-03
								6.6656320041650000E-03]';
                case 73 % error is of order O(h^20)
                    pt = [3.3333333333333300E-01	3.3333333333333300E-01
                        2.0780025853987000E-02	4.8960998707300600E-01
                        4.8960998707300600E-01	4.8960998707300600E-01
                        4.8960998707300600E-01	2.0780025853987000E-02
                        9.0926214604215000E-02	4.5453689269789200E-01
                        4.5453689269789200E-01	4.5453689269789200E-01
                        4.5453689269789200E-01	9.0926214604215000E-02
                        1.9716663870113700E-01	4.0141668064943000E-01
                        4.0141668064943000E-01	4.0141668064943000E-01
                        4.0141668064943000E-01	1.9716663870113700E-01
                        4.8889669119380500E-01	2.5555165440309800E-01
                        2.5555165440309800E-01	2.5555165440309800E-01
                        2.5555165440309800E-01	4.8889669119380500E-01
                        6.4584411569574000E-01	1.7707794215212900E-01
                        1.7707794215212900E-01	1.7707794215212900E-01
                        1.7707794215212900E-01	6.4584411569574000E-01
                        7.7987789354409500E-01	1.1006105322795200E-01
                        1.1006105322795200E-01	1.1006105322795200E-01
                        1.1006105322795200E-01	7.7987789354409500E-01
                        8.8894275149632000E-01	5.5528624251839900E-02
                        5.5528624251839900E-02	5.5528624251839900E-02
                        5.5528624251839900E-02	8.8894275149632000E-01
                        9.7475627244554200E-01	1.2621863777229000E-02
                        1.2621863777229000E-02	1.2621863777229000E-02
                        1.2621863777229000E-02	9.7475627244554200E-01
                        3.6114178484120000E-03	3.9575478735694300E-01
                        3.9575478735694300E-01	6.0063379479464400E-01
                        6.0063379479464400E-01	3.6114178484120000E-03
                        3.9575478735694300E-01	3.6114178484120000E-03
                        6.0063379479464400E-01	3.9575478735694300E-01
                        3.6114178484120000E-03	6.0063379479464400E-01
                        1.3446675453078000E-01	3.0792998388043600E-01
                        3.0792998388043600E-01	5.5760326158878400E-01
                        5.5760326158878400E-01	1.3446675453078000E-01
                        3.0792998388043600E-01	1.3446675453078000E-01
                        5.5760326158878400E-01	3.0792998388043600E-01
                        1.3446675453078000E-01	5.5760326158878400E-01
                        1.4446025776115000E-02	2.6456694840652000E-01
                        2.6456694840652000E-01	7.2098702581736400E-01
                        7.2098702581736400E-01	1.4446025776115000E-02
                        2.6456694840652000E-01	1.4446025776115000E-02
                        7.2098702581736400E-01	2.6456694840652000E-01
                        1.4446025776115000E-02	7.2098702581736400E-01
                        4.6933578838178000E-02	3.5853935220595100E-01
                        3.5853935220595100E-01	5.9452706895587100E-01
                        5.9452706895587100E-01	4.6933578838178000E-02
                        3.5853935220595100E-01	4.6933578838178000E-02
                        5.9452706895587100E-01	3.5853935220595100E-01
                        4.6933578838178000E-02	5.9452706895587100E-01
                        2.8611203505669900E-03	1.5780740596859400E-01
                        1.5780740596859400E-01	8.3933147368083900E-01
                        8.3933147368083900E-01	2.8611203505669900E-03
                        1.5780740596859400E-01	2.8611203505669900E-03
                        8.3933147368083900E-01	1.5780740596859400E-01
                        2.8611203505669900E-03	8.3933147368083900E-01
                        2.2386142409791600E-01	7.5050596975911000E-02
                        7.5050596975911000E-02	7.0108797892617300E-01
                        7.0108797892617300E-01	2.2386142409791600E-01
                        7.5050596975911000E-02	2.2386142409791600E-01
                        7.0108797892617300E-01	7.5050596975911000E-02
                        2.2386142409791600E-01	7.0108797892617300E-01
                        3.4647074816760000E-02	1.4242160111338300E-01
                        1.4242160111338300E-01	8.2293132406985700E-01
                        8.2293132406985700E-01	3.4647074816760000E-02
                        1.4242160111338300E-01	3.4647074816760000E-02
                        8.2293132406985700E-01	1.4242160111338300E-01
                        3.4647074816760000E-02	8.2293132406985700E-01
                        1.0161119296278000E-02	6.5494628082938000E-02
                        6.5494628082938000E-02	9.2434425262078400E-01
                        9.2434425262078400E-01	1.0161119296278000E-02
                        6.5494628082938000E-02	1.0161119296278000E-02
                        9.2434425262078400E-01	6.5494628082938000E-02
                        1.0161119296278000E-02	9.2434425262078400E-01]';
                weight = 0.5*[3.2906331388919000E-02
                            1.0330731891272000E-02
                            1.0330731891272000E-02
                            1.0330731891272000E-02
                            2.2387247263016000E-02
                            2.2387247263016000E-02
                            2.2387247263016000E-02
                            3.0266125869467900E-02
                            3.0266125869467900E-02
                            3.0266125869467900E-02
                            3.0490967802198000E-02
                            3.0490967802198000E-02
                            3.0490967802198000E-02
                            2.4159212741641000E-02
                            2.4159212741641000E-02
                            2.4159212741641000E-02
                            1.6050803586800900E-02
                            1.6050803586800900E-02
                            1.6050803586800900E-02
                            8.0845802617839900E-03
                            8.0845802617839900E-03
                            8.0845802617839900E-03
                            2.0793620274849900E-03
                            2.0793620274849900E-03
                            2.0793620274849900E-03
                            3.8848769049810000E-03
                            3.8848769049810000E-03
                            3.8848769049810000E-03
                            3.8848769049810000E-03
                            3.8848769049810000E-03
                            3.8848769049810000E-03
                            2.5574160612022000E-02
                            2.5574160612022000E-02
                            2.5574160612022000E-02
                            2.5574160612022000E-02
                            2.5574160612022000E-02
                            2.5574160612022000E-02
                            8.8809035733379900E-03
                            8.8809035733379900E-03
                            8.8809035733379900E-03
                            8.8809035733379900E-03
                            8.8809035733379900E-03
                            8.8809035733379900E-03
                            1.6124546761731000E-02
                            1.6124546761731000E-02
                            1.6124546761731000E-02
                            1.6124546761731000E-02
                            1.6124546761731000E-02
                            1.6124546761731000E-02
                            2.4919418174909900E-03
                            2.4919418174909900E-03
                            2.4919418174909900E-03
                            2.4919418174909900E-03
                            2.4919418174909900E-03
                            2.4919418174909900E-03
                            1.8242840118951000E-02
                            1.8242840118951000E-02
                            1.8242840118951000E-02
                            1.8242840118951000E-02
                            1.8242840118951000E-02
                            1.8242840118951000E-02
                            1.0258563736198900E-02
                            1.0258563736198900E-02
                            1.0258563736198900E-02
                            1.0258563736198900E-02
                            1.0258563736198900E-02
                            1.0258563736198900E-02
                            3.7999288553020000E-03
                            3.7999288553020000E-03
                            3.7999288553020000E-03
                            3.7999288553020000E-03
                            3.7999288553020000E-03
                            3.7999288553020000E-03]';


                otherwise
                    error('gauss_points_and_weights: number of gauss points unavailable')
            end
        case 3
            switch npt
                case 4 % error is of order O(h^3)
                    r1 = 0.585410196624969;
                    r2 = 0.138196601125011;
                    w1 = 1./24.;
                    pt = [r1,r2,r2;
                          r2,r1,r2;
                          r2,r2,r1;
                          r2,r2,r2]';
                    weight = [w1,w1,w1,w1];
                case 5 % error is of order O(h^4)
                    r1 = 0.5;
                    r2 = 1./6.;
                    w1 = -0.8/6.;
                    w2 = 0.45/6.;
                    pt = [0.25,0.25,0.25;
                            r1,r2,r2;
                            r2,r1,r2;
                            r2,r2,r1;
                            r2,r2,r2]';
                    weight = [w1,w2,w2,w2,w2];
                case 11 % error is of order O(h^5)
                    w1 = -0.013155555555556;
                    w2 = 0.007622222222222;
                    w3 = 0.024888888888889;
                    r1 = 0.785714285714286;
                    r2 = 0.071428571428571;
                    r3 = 0.399403576166799;
                    r4 = 0.100596423833201;
                    pt = [0.25,0.25,0.25;
                                r1,r2,r2;
                                r2,r1,r2;
                                r2,r2,r1;
                                r2,r2,r2;
                                r3,r4,r4;
                                r4,r3,r4;
                                r4,r4,r3;
                                r4,r3,r3;
                                r3,r4,r3;
                                r3,r3,r4]';
                    weight = [w1,w2,w2,w2,w2,w3,w3,w3,w3,w3,w3];   
                case 15 % error is of order O(h^6)
                    w1 = 0.030283678097089;
                    w2 = 0.006026785714286;
                    w3 = 0.011645249086029;
                    w4 = 0.010949141561386;
                    r1 = 1./3.;
                    r2 = 0.727272727272727;
                    r3 = 0.090909090909091;
                    r4 = 0.066550153573664;
                    r5 = 0.433449846426336;
                    pt = [0.25,0.25,0.25;
                                0,r1,r1;
                                r1,0,r1;
                                r1,r1,0;
                                r1,r1,r1;
                                r2,r3,r3;
                                r3,r2,r3;
                                r3,r3,r2;
                                r3,r3,r3;
                                r4,r5,r5;
                                r5,r4,r5;
                                r5,r5,r4;
                                r5,r4,r4;
                                r4,r5,r4;
                                r4,r4,r5]';
                    weight = [w1,w2,w2,w2,w2,w3,w3,w3,w3,w4,w4,w4,w4,w4,w4];



                otherwise
                    error('gauss_points_and_weights: number of gauss points unavailable')    
            end  
        otherwise
            error('gauss_points_and_weights: dimension unavailable')
    end
else     
    [pt1,weight1]=lgwt(npt,-1,1);
    switch dimension
        case 1
            pt = pt1';
            weight = weight1';
        case 2
            pt= zeros(2,npt*npt);
            weight = zeros(1,npt*npt);
            for i=1:npt
                for j=1:npt
                    k=npt*(i-1)+j;
                    pt(1:2,k)=[pt1(i);pt1(j)];
                    weight(k)=weight1(i)*weight1(j);
                end
            end                 
        case 3
            nptsq = npt*npt;
            pt = zeros(3,npt*nptsq);
            weight = zeros(1,npt*nptsq);
            for i=1:npt
                for j=1:npt
                    for k=1:npt
                        p = nptsq*(i-1)+npt*(j-1)+k;
                        pt(1:3,p) = [pt1(i);pt1(j);pt1(k)];
                        weight(p) = weight1(i)*weight1(j)*weight1(k);
                    end
                end
            end
        otherwise
            error('gauss_points_and_weights: number of gauss points unavailable')
     end
end  
if (formatID == 0)
    gauss.pt = pt;
    gauss.weight = weight;
    gauss.npt = numel(gauss.weight);
elseif (formatID == 1)    
    gauss = [pt;weight];
end

end
