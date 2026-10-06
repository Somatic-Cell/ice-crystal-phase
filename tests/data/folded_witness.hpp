#pragma once
// Generated from the supplied GPU snapshot, not a synthetic optical result.
// Source SHA256: b565245e04fc0576f9ebeaea843cf623ab0f6ec38e61c4ea5efbb12226c8f8f3
#include <rainbow/patch_optical_evaluator.hpp>
#include <bit>
#include <vector>
namespace folded_tests {
using namespace rainbow;
inline float F(std::uint32_t x) { return std::bit_cast<float>(x); }
struct Witness {
 std::vector<OutgoingVertex> vertices;
 std::vector<OutgoingPatch> patches;
 std::vector<std::uint32_t> used_patches;
 std::vector<Vec3> directions;
 std::vector<PatchQuerySummary> summaries;
 std::vector<PatchQueryHit> hits;
 std::vector<std::uint64_t> offsets;
 std::vector<std::uint32_t> original_direction_ids;
 FocalPhaseConfig focal{};
 Vec3 incident{},basis{};
 PatchOpticsParams view() const {
  PatchOpticsParams p{};p.vertices=vertices.data();p.vertex_count=static_cast<std::uint32_t>(vertices.size());
  p.patches=patches.data();p.patch_count=static_cast<std::uint32_t>(patches.size());
  p.directions=directions.data();p.direction_count=static_cast<std::uint32_t>(directions.size());
  p.hits=hits.data();p.hit_storage_count=hits.size();p.offsets=offsets.data();p.summaries=summaries.data();
  p.incident_direction=incident;p.incident_basis_x=basis;return p;
 }
};
inline Witness make_witness() {
 Witness a;
 a.vertices.resize(66564);a.patches.resize(43622);
 a.incident={F(1064341426u),F(3199147331u),F(0u)};a.basis={F(0u),F(0u),F(1065353216u)};
 a.focal={129,129,F(1065848027u),0,{0,0,0,0}};
 a.vertices[5852]={{F(3208481698u),F(1058710294u),F(3197109673u)},{F(3199065421u),F(1062086365u),F(3203978960u)},{F(1045672133u),F(1058105179u),F(1061941641u)},{{F(3188253044u),F(0u)},{F(1032845773u),F(0u)}},{3452,F(1064770014u)},VertexStatus::Valid,0};
 a.vertices[5853]={{F(3208558668u),F(1058738309u),F(3196554459u)},{F(3199556842u),F(1062237822u),F(3203099782u)},{F(1045135446u),F(1057736545u),F(1062230403u)},{{F(3188378368u),F(0u)},{F(1032579982u),F(0u)}},{3439,F(1062728759u)},VertexStatus::Valid,0};
 a.vertices[5981]={{F(3208678676u),F(1058486565u),F(3197109673u)},{F(3200075632u),F(1061824047u),F(3204102203u)},{F(1046087628u),F(1058390569u),F(1061705681u)},{{F(3188005258u),F(0u)},{F(1033272694u),F(0u)}},{3432,F(1060887554u)},VertexStatus::Valid,0};
 a.vertices[5982]={{F(3208755129u),F(1058514392u),F(3196554459u)},{F(3200569612u),F(1061972837u),F(3203215359u)},{F(1045553131u),F(1058023440u),F(1062007194u)},{{F(3188139530u),F(0u)},{F(1033003987u),F(0u)}},{3419,F(1059990052u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3315];p.patch_id=5807;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129619244u);p.vertex_indices[0]=5852;p.vertex_indices[1]=5853;p.vertex_indices[2]=5981;p.vertex_indices[3]=5982;a.used_patches.push_back(3315);}
 a.vertices[5885]={{F(3208558668u),F(1058738309u),F(1049070811u)},{F(3199556842u),F(1062237822u),F(1055616134u)},{F(3192619094u),F(3205220193u),F(1062230403u)},{{F(3188378368u),F(0u)},{F(3180063630u),F(0u)}},{3439,F(1062728759u)},VertexStatus::Valid,0};
 a.vertices[5886]={{F(3208481698u),F(1058710294u),F(1049626025u)},{F(3199065421u),F(1062086365u),F(1056495312u)},{F(3193155781u),F(3205588827u),F(1061941641u)},{{F(3188253044u),F(0u)},{F(3180329421u),F(0u)}},{3452,F(1064770014u)},VertexStatus::Valid,0};
 a.vertices[6014]={{F(3208755129u),F(1058514392u),F(1049070811u)},{F(3200569612u),F(1061972837u),F(1055731711u)},{F(3193036779u),F(3205507088u),F(1062007194u)},{{F(3188139530u),F(0u)},{F(3180487635u),F(0u)}},{3419,F(1059990052u)},VertexStatus::Valid,0};
 a.vertices[6015]={{F(3208678676u),F(1058486565u),F(1049626025u)},{F(3200075632u),F(1061824047u),F(1056618555u)},{F(3193571276u),F(3205874217u),F(1061705681u)},{{F(3188005258u),F(0u)},{F(3180756342u),F(0u)}},{3432,F(1060887554u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3348];p.patch_id=5840;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129619244u);p.vertex_indices[0]=5885;p.vertex_indices[1]=5886;p.vertex_indices[2]=6014;p.vertex_indices[3]=6015;a.used_patches.push_back(3348);}
 a.vertices[5983]={{F(3208826587u),F(1058540400u),F(3195938843u)},{F(3201033734u),F(1062112360u),F(3202313654u)},{F(1044981643u),F(1057630903u),F(1062309971u)},{{F(3188276703u),F(0u)},{F(1032694199u),F(0u)}},{3407,F(1056431186u)},VertexStatus::Valid,0};
 a.vertices[6111]={{F(3208946164u),F(1058288500u),F(3196554459u)},{F(3201556020u),F(1061691481u),F(3203324915u)},{F(1045989516u),F(1058323179u),F(1061762410u)},{{F(3187895427u),F(0u)},{F(1033440120u),F(0u)}},{3400,F(1055513846u)},VertexStatus::Valid,0};
 a.vertices[6112]={{F(3209017179u),F(1058314348u),F(3195938843u)},{F(3202022485u),F(1061828559u),F(3202416028u)},{F(1045419911u),F(1057931935u),F(1062079536u)},{{F(3188042255u),F(0u)},{F(1033127539u),F(0u)}},{3388,F(1050845626u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3432];p.patch_id=5936;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129623992u);p.vertex_indices[0]=5982;p.vertex_indices[1]=5983;p.vertex_indices[2]=6111;p.vertex_indices[3]=6112;a.used_patches.push_back(3432);}
 a.vertices[6013]={{F(3208826587u),F(1058540400u),F(1048455195u)},{F(3201033734u),F(1062112360u),F(1054830006u)},{F(3192465291u),F(3205114551u),F(1062309971u)},{{F(3188276703u),F(0u)},{F(3180177847u),F(0u)}},{3407,F(1056431186u)},VertexStatus::Valid,0};
 a.vertices[6142]={{F(3209017179u),F(1058314348u),F(1048455195u)},{F(3202022485u),F(1061828559u),F(1054932380u)},{F(3192903559u),F(3205415583u),F(1062079536u)},{{F(3188042255u),F(0u)},{F(3180611187u),F(0u)}},{3388,F(1050845626u)},VertexStatus::Valid,0};
 a.vertices[6143]={{F(3208946164u),F(1058288500u),F(1049070811u)},{F(3201556020u),F(1061691481u),F(1055841267u)},{F(3193473164u),F(3205806827u),F(1061762410u)},{{F(3187895427u),F(0u)},{F(3180923768u),F(0u)}},{3400,F(1055513846u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3463];p.patch_id=5967;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129623992u);p.vertex_indices[0]=6013;p.vertex_indices[1]=6014;p.vertex_indices[2]=6142;p.vertex_indices[3]=6143;a.used_patches.push_back(3463);}
 a.vertices[6107]={{F(3208611709u),F(1058166769u),F(3198775314u)},{F(3199390212u),F(1061051738u),F(3205625014u)},{F(1047898679u),F(1059634526u),F(1060533574u)},{{F(3187031477u),F(0u)},{F(1034273635u),F(0u)}},{3457,F(1059990900u)},VertexStatus::Valid,0};
 a.vertices[6108]={{F(3208703051u),F(1058200014u),F(3198220100u)},{F(3199976615u),F(1061225499u),F(3205203411u)},{F(1047472249u),F(1059341624u),F(1060832538u)},{{F(3187287877u),F(0u)},{F(1034125460u),F(0u)}},{3442,F(1026709840u)},VertexStatus::Valid,0};
 a.vertices[6236]={{F(3208799436u),F(1057939673u),F(3198775314u)},{F(3200339034u),F(1060766000u),F(3205690815u)},{F(1048319954u),F(1059923886u),F(1060222508u)},{{F(3186466490u),F(0u)},{F(1034710703u),F(0u)}},{3439,F(1039533264u)},VertexStatus::Valid,0};
 a.vertices[6237]={{F(3208890222u),F(1057972716u),F(3198220100u)},{F(3200928363u),F(1060936638u),F(3205265651u)},{F(1047903970u),F(1059638159u),F(1060529766u)},{{F(3186736462u),F(0u)},{F(1034566869u),F(0u)}},{3423,F(1057903170u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3544];p.patch_id=6060;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129631646u);p.vertex_indices[0]=6107;p.vertex_indices[1]=6108;p.vertex_indices[2]=6236;p.vertex_indices[3]=6237;a.used_patches.push_back(3544);}
 a.vertices[6113]={{F(3209083301u),F(1058338414u),F(3194828415u)},{F(3202458916u),F(1061956593u),F(3201493398u)},{F(1044808987u),F(1057512310u),F(1062397659u)},{{F(3188191661u),F(0u)},{F(1032770982u),F(0u)}},{3377,F(1016629376u)},VertexStatus::Valid,0};
 a.vertices[6241]={{F(3209202478u),F(1058086368u),F(3195938843u)},{F(3202984340u),F(1061529015u),F(3202512887u)},{F(1045879327u),F(1058247494u),F(1061825370u)},{{F(3187800889u),F(0u)},{F(1033575457u),F(0u)}},{3370,F(1032987080u)},VertexStatus::Valid,0};
 a.vertices[6242]={{F(3209268219u),F(1058110296u),F(3194828415u)},{F(3203422914u),F(1061654809u),F(3201583531u)},{F(1045269815u),F(1057828839u),F(1062159746u)},{{F(3187960851u),F(0u)},{F(1033216037u),F(0u)}},{3358,F(1062341285u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3549];p.patch_id=6065;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129628594u);p.vertex_indices[0]=6112;p.vertex_indices[1]=6113;p.vertex_indices[2]=6241;p.vertex_indices[3]=6242;a.used_patches.push_back(3549);}
 a.vertices[6141]={{F(3209083301u),F(1058338414u),F(1047344767u)},{F(3202458916u),F(1061956593u),F(1054009750u)},{F(3192292635u),F(3204995958u),F(1062397659u)},{{F(3188191661u),F(0u)},{F(3180254630u),F(0u)}},{3377,F(1016629376u)},VertexStatus::Valid,0};
 a.vertices[6270]={{F(3209268219u),F(1058110296u),F(1047344767u)},{F(3203422914u),F(1061654809u),F(1054099883u)},{F(3192753463u),F(3205312487u),F(1062159746u)},{{F(3187960851u),F(0u)},{F(3180699685u),F(0u)}},{3358,F(1062341285u)},VertexStatus::Valid,0};
 a.vertices[6271]={{F(3209202478u),F(1058086368u),F(1048455195u)},{F(3202984340u),F(1061529015u),F(1055029239u)},{F(3193362975u),F(3205731142u),F(1061825370u)},{{F(3187800889u),F(0u)},{F(3181059105u),F(0u)}},{3370,F(1032987080u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3578];p.patch_id=6094;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129628594u);p.vertex_indices[0]=6141;p.vertex_indices[1]=6142;p.vertex_indices[2]=6270;p.vertex_indices[3]=6271;a.used_patches.push_back(3578);}
 a.vertices[6146]={{F(3208703051u),F(1058200014u),F(1050736452u)},{F(3199976615u),F(1061225499u),F(1057719763u)},{F(3194955897u),F(3206825272u),F(1060832538u)},{{F(3187287877u),F(0u)},{F(3181609108u),F(0u)}},{3442,F(1026709840u)},VertexStatus::Valid,0};
 a.vertices[6147]={{F(3208611709u),F(1058166769u),F(1051291666u)},{F(3199390212u),F(1061051738u),F(1058141366u)},{F(3195382327u),F(3207118174u),F(1060533574u)},{{F(3187031477u),F(0u)},{F(3181757283u),F(0u)}},{3457,F(1059990900u)},VertexStatus::Valid,0};
 a.vertices[6275]={{F(3208890222u),F(1057972716u),F(1050736452u)},{F(3200928363u),F(1060936638u),F(1057782003u)},{F(3195387618u),F(3207121807u),F(1060529766u)},{{F(3186736462u),F(0u)},{F(3182050517u),F(0u)}},{3423,F(1057903170u)},VertexStatus::Valid,0};
 a.vertices[6276]={{F(3208799436u),F(1057939673u),F(1051291666u)},{F(3200339034u),F(1060766000u),F(1058207167u)},{F(3195803602u),F(3207407534u),F(1060222508u)},{{F(3186466490u),F(0u)},{F(3182194351u),F(0u)}},{3439,F(1039533264u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3583];p.patch_id=6099;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129631646u);p.vertex_indices[0]=6146;p.vertex_indices[1]=6147;p.vertex_indices[2]=6275;p.vertex_indices[3]=6276;a.used_patches.push_back(3583);}
 a.vertices[6238]={{F(3208975834u),F(1058003876u),F(3197664886u)},{F(3201487578u),F(1061098206u),F(3204831487u)},{F(1047454532u),F(1059329454u),F(1060844631u)},{{F(3187019014u),F(0u)},{F(1034382891u),F(0u)}},{3408,F(1063417106u)},VertexStatus::Valid,0};
 a.vertices[6366]={{F(3209072065u),F(1057743479u),F(3198220100u)},{F(3201852959u),F(1060632443u),F(3205324391u)},{F(1048347607u),F(1059942880u),F(1060201513u)},{{F(3186168818u),F(0u)},{F(1035014010u),F(0u)}},{3405,F(1065272731u)},VertexStatus::Valid,0};
 a.vertices[6367]={{F(3209157190u),F(1057774462u),F(3197664886u)},{F(3202414872u),F(1060791122u),F(3204886915u)},{F(1047909783u),F(1059642152u),F(1060525581u)},{{F(3186465566u),F(0u)},{F(1034835787u),F(0u)}},{3391,F(1053149118u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3662];p.patch_id=6189;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129635673u);p.vertex_indices[0]=6237;p.vertex_indices[1]=6238;p.vertex_indices[2]=6366;p.vertex_indices[3]=6367;a.used_patches.push_back(3662);}
 a.vertices[6274]={{F(3208975834u),F(1058003876u),F(1050181238u)},{F(3201487578u),F(1061098206u),F(1057347839u)},{F(3194938180u),F(3206813102u),F(1060844631u)},{{F(3187019014u),F(0u)},{F(3181866539u),F(0u)}},{3408,F(1063417106u)},VertexStatus::Valid,0};
 a.vertices[6403]={{F(3209157190u),F(1057774462u),F(1050181238u)},{F(3202414872u),F(1060791122u),F(1057403267u)},{F(3195393431u),F(3207125800u),F(1060525581u)},{{F(3186465566u),F(0u)},{F(3182319435u),F(0u)}},{3391,F(1053149118u)},VertexStatus::Valid,0};
 a.vertices[6404]={{F(3209072065u),F(1057743479u),F(1050736452u)},{F(3201852959u),F(1060632443u),F(1057840743u)},{F(3195831255u),F(3207426528u),F(1060201513u)},{{F(3186168818u),F(0u)},{F(3182497658u),F(0u)}},{3405,F(1065272731u)},VertexStatus::Valid,0};
 {auto& p=a.patches[3699];p.patch_id=6226;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129635674u);p.vertex_indices[0]=6274;p.vertex_indices[1]=6275;p.vertex_indices[2]=6403;p.vertex_indices[3]=6404;a.used_patches.push_back(3699);}
 a.vertices[7003]={{F(3209050369u),F(1055552325u),F(3202661810u)},{F(3200342044u),F(1057448892u),F(3208811441u)},{F(1050517660u),F(1062767089u),F(1054803736u)},{{F(3180968783u),F(0u)},{F(1036905500u),F(0u)}},{3479,F(1056500402u)},VertexStatus::Valid,0};
 a.vertices[7004]={{F(3209176525u),F(1055644159u),F(3202106596u)},{F(3201164253u),F(1057659238u),F(3208447266u)},{F(1050440762u),F(1062661451u),F(1055259502u)},{{F(3181154151u),F(0u)},{F(1037065590u),F(0u)}},{3457,F(1058957848u)},VertexStatus::Valid,0};
 a.vertices[7132]={{F(3209205213u),F(1055074196u),F(3202661810u)},{F(3201070624u),F(1057091262u),F(3208868945u)},{F(1050669316u),F(1062975425u),F(1053849580u)},{{F(3180254666u),F(0u)},{F(1037232395u),F(0u)}},{3466,F(1063065749u)},VertexStatus::Valid,0};
 a.vertices[7133]={{F(3209330871u),F(1055165667u),F(3202106596u)},{F(3201896427u),F(1057297524u),F(3208502300u)},{F(1050600987u),F(1062881559u),F(1054289196u)},{{F(3180436776u),F(0u)},{F(1037403683u),F(0u)}},{3445,F(1029884976u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4369];p.patch_id=6949;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129676716u);p.vertex_indices[0]=7003;p.vertex_indices[1]=7004;p.vertex_indices[2]=7132;p.vertex_indices[3]=7133;a.used_patches.push_back(4369);}
 a.vertices[7005]={{F(3209296879u),F(1055731769u),F(3201551382u)},{F(3201955591u),F(1057861209u),F(3208069967u)},{F(1050356472u),F(1062545658u),F(1055740020u)},{{F(3181356565u),F(0u)},{F(1037195503u),F(0u)}},{3436,F(1061436440u)},VertexStatus::Valid,0};
 a.vertices[7134]={{F(3209450755u),F(1055252935u),F(3201551382u)},{F(3202691200u),F(1057495603u),F(3208122581u)},{F(1050525871u),F(1062778370u),F(1054754015u)},{{F(3180635878u),F(0u)},{F(1037545981u),F(0u)}},{3424,F(1048854137u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4370];p.patch_id=6950;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129675254u);p.vertex_indices[0]=7004;p.vertex_indices[1]=7005;p.vertex_indices[2]=7133;p.vertex_indices[3]=7134;a.used_patches.push_back(4370);}
 a.vertices[7006]={{F(3209411572u),F(1055815259u),F(3200996168u)},{F(3202716031u),F(1058054867u),F(3207680240u)},{F(1050263882u),F(1062418465u),F(1056246938u)},{{F(3181576800u),F(0u)},{F(1037294837u),F(0u)}},{3416,F(1063733532u)},VertexStatus::Valid,0};
 a.vertices[7135]={{F(3209565008u),F(1055336105u),F(3200996168u)},{F(3203454921u),F(1057685555u),F(3207730478u)},{F(1050443098u),F(1062664660u),F(1055245912u)},{{F(3180852879u),F(0u)},{F(1037658940u),F(0u)}},{3404,F(1055163560u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4371];p.patch_id=6951;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129673877u);p.vertex_indices[0]=7005;p.vertex_indices[1]=7006;p.vertex_indices[2]=7134;p.vertex_indices[3]=7135;a.used_patches.push_back(4371);}
 a.vertices[7007]={{F(3209520738u),F(1055894725u),F(3200440955u)},{F(3203445538u),F(1058240265u),F(3207278752u)},{F(1050161956u),F(1062278446u),F(1056781970u)},{{F(3181815682u),F(0u)},{F(1037362858u),F(0u)}},{3398,F(1002179104u)},VertexStatus::Valid,0};
 a.vertices[7136]={{F(3209673766u),F(1055415274u),F(3200440955u)},{F(3204187578u),F(1057867428u),F(3207326650u)},{F(1050351661u),F(1062539051u),F(1055766877u)},{{F(3181088777u),F(0u)},{F(1037741858u),F(0u)}},{3385,F(1058555873u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4372];p.patch_id=6952;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129672617u);p.vertex_indices[0]=7006;p.vertex_indices[1]=7007;p.vertex_indices[2]=7135;p.vertex_indices[3]=7136;a.used_patches.push_back(4372);}
 a.vertices[7053]={{F(3209520738u),F(1055894725u),F(1052957307u)},{F(3203445538u),F(1058240265u),F(1059795104u)},{F(3197645604u),F(3209762094u),F(1056781970u)},{{F(3181815682u),F(0u)},{F(3184846506u),F(0u)}},{3398,F(1002179104u)},VertexStatus::Valid,0};
 a.vertices[7054]={{F(3209411572u),F(1055815259u),F(1053512520u)},{F(3202716031u),F(1058054867u),F(1060196592u)},{F(3197747530u),F(3209902113u),F(1056246938u)},{{F(3181576800u),F(0u)},{F(3184778485u),F(0u)}},{3416,F(1063733532u)},VertexStatus::Valid,0};
 a.vertices[7182]={{F(3209673766u),F(1055415274u),F(1052957307u)},{F(3204187578u),F(1057867428u),F(1059843002u)},{F(3197835309u),F(3210022699u),F(1055766877u)},{{F(3181088777u),F(0u)},{F(3185225506u),F(0u)}},{3385,F(1058555873u)},VertexStatus::Valid,0};
 a.vertices[7183]={{F(3209565008u),F(1055336105u),F(1053512520u)},{F(3203454921u),F(1057685555u),F(1060246830u)},{F(3197926746u),F(3210148308u),F(1055245912u)},{{F(3180852879u),F(0u)},{F(3185142588u),F(0u)}},{3404,F(1055163560u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4419];p.patch_id=6999;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129672616u);p.vertex_indices[0]=7053;p.vertex_indices[1]=7054;p.vertex_indices[2]=7182;p.vertex_indices[3]=7183;a.used_patches.push_back(4419);}
 a.vertices[7055]={{F(3209296879u),F(1055731769u),F(1054067734u)},{F(3201955591u),F(1057861209u),F(1060586319u)},{F(3197840120u),F(3210029306u),F(1055740020u)},{{F(3181356565u),F(0u)},{F(3184679151u),F(0u)}},{3436,F(1061436440u)},VertexStatus::Valid,0};
 a.vertices[7184]={{F(3209450755u),F(1055252935u),F(1054067734u)},{F(3202691200u),F(1057495603u),F(1060638933u)},{F(3198009519u),F(3210262018u),F(1054754015u)},{{F(3180635878u),F(0u)},{F(3185029629u),F(0u)}},{3424,F(1048854137u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4420];p.patch_id=7000;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129673877u);p.vertex_indices[0]=7054;p.vertex_indices[1]=7055;p.vertex_indices[2]=7183;p.vertex_indices[3]=7184;a.used_patches.push_back(4420);}
 a.vertices[7056]={{F(3209176525u),F(1055644159u),F(1054622948u)},{F(3201164253u),F(1057659238u),F(1060963618u)},{F(3197924410u),F(3210145099u),F(1055259502u)},{{F(3181154151u),F(0u)},{F(3184549238u),F(0u)}},{3457,F(1058957848u)},VertexStatus::Valid,0};
 a.vertices[7185]={{F(3209330871u),F(1055165667u),F(1054622948u)},{F(3201896427u),F(1057297524u),F(1061018652u)},{F(3198084635u),F(3210365207u),F(1054289196u)},{{F(3180436776u),F(0u)},{F(3184887331u),F(0u)}},{3445,F(1029884976u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4421];p.patch_id=7001;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129675254u);p.vertex_indices[0]=7055;p.vertex_indices[1]=7056;p.vertex_indices[2]=7184;p.vertex_indices[3]=7185;a.used_patches.push_back(4421);}
 a.vertices[7057]={{F(3209050369u),F(1055552325u),F(1055178162u)},{F(3200342044u),F(1057448892u),F(1061327793u)},{F(3198001308u),F(3210250737u),F(1054803736u)},{{F(3180968783u),F(0u)},{F(3184389148u),F(0u)}},{3479,F(1056500402u)},VertexStatus::Valid,0};
 a.vertices[7186]={{F(3209205213u),F(1055074196u),F(1055178162u)},{F(3201070624u),F(1057091262u),F(1061385297u)},{F(3198152964u),F(3210459073u),F(1053849580u)},{{F(3180254666u),F(0u)},{F(3184716043u),F(0u)}},{3466,F(1063065749u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4422];p.patch_id=7002;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129676716u);p.vertex_indices[0]=7056;p.vertex_indices[1]=7057;p.vertex_indices[2]=7185;p.vertex_indices[3]=7186;a.used_patches.push_back(4422);}
 a.vertices[7262]={{F(3209480138u),F(1054683480u),F(3202106596u)},{F(3202599020u),F(1056882787u),F(3208552839u)},{F(1050754903u),F(1063092999u),F(1053273820u)},{{F(3179701788u),F(0u)},{F(1037725285u),F(0u)}},{3433,F(1053239810u)},VertexStatus::Valid,0};
 a.vertices[7263]={{F(3209599602u),F(1054770442u),F(3201551382u)},{F(3203397188u),F(1057117930u),F(3208170903u)},{F(1050689093u),F(1063002594u),F(1053719121u)},{{F(3179895368u),F(0u)},{F(1037880517u),F(0u)}},{3412,F(1059291986u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4490];p.patch_id=7078;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129683344u);p.vertex_indices[0]=7133;p.vertex_indices[1]=7134;p.vertex_indices[2]=7262;p.vertex_indices[3]=7263;a.used_patches.push_back(4490);}
 a.vertices[7264]={{F(3209713458u),F(1054853323u),F(3200996168u)},{F(3204164146u),F(1057304220u),F(3207776623u)},{F(1050616347u),F(1062902660u),F(1054191822u)},{{F(3180106751u),F(0u)},{F(1038007750u),F(0u)}},{3392,F(1063183919u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4491];p.patch_id=7079;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129681826u);p.vertex_indices[0]=7134;p.vertex_indices[1]=7135;p.vertex_indices[2]=7263;p.vertex_indices[3]=7264;a.used_patches.push_back(4491);}
 a.vertices[7312]={{F(3209713458u),F(1054853323u),F(1053512520u)},{F(3204164146u),F(1057304220u),F(1060292975u)},{F(3198099995u),F(3210386308u),F(1054191822u)},{{F(3180106751u),F(0u)},{F(3185491398u),F(0u)}},{3392,F(1063183919u)},VertexStatus::Valid,0};
 a.vertices[7313]={{F(3209599602u),F(1054770442u),F(1054067734u)},{F(3203397188u),F(1057117930u),F(1060687255u)},{F(3198172741u),F(3210486242u),F(1053719121u)},{{F(3179895368u),F(0u)},{F(3185364165u),F(0u)}},{3412,F(1059291986u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4540];p.patch_id=7128;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129681826u);p.vertex_indices[0]=7183;p.vertex_indices[1]=7184;p.vertex_indices[2]=7312;p.vertex_indices[3]=7313;a.used_patches.push_back(4540);}
 a.vertices[7314]={{F(3209480138u),F(1054683480u),F(1054622948u)},{F(3202599020u),F(1056882787u),F(1061069191u)},{F(3198238551u),F(3210576647u),F(1053273820u)},{{F(3179701788u),F(0u)},{F(3185208933u),F(0u)}},{3433,F(1053239810u)},VertexStatus::Valid,0};
 {auto& p=a.patches[4541];p.patch_id=7129;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129683344u);p.vertex_indices[0]=7184;p.vertex_indices[1]=7185;p.vertex_indices[2]=7313;p.vertex_indices[3]=7314;a.used_patches.push_back(4541);}
 a.vertices[10364]={{F(3212235016u),F(1032694847u),F(3198775314u)},{F(3209293391u),F(3190062923u),F(3206013509u)},{F(1049548714u),F(1061436010u),F(3205764623u)},{{F(1035924910u),F(0u)},{F(1036876244u),F(0u)}},{3312,F(1051511494u)},VertexStatus::Valid,0};
 a.vertices[10365]={{F(3212324233u),F(1032954624u),F(3198220100u)},{F(3209629709u),F(3189731826u),F(3205573011u)},{F(1049386488u),F(1061213154u),F(3206090674u)},{{F(1036232646u),F(0u)},{F(1036795618u),F(0u)}},{3296,F(1048844928u)},VertexStatus::Valid,0};
 a.vertices[10493]={{F(3212262181u),F(1029022329u),F(3198775314u)},{F(3209199734u),F(3192368336u),F(3205964727u)},{F(1049306391u),F(1061103122u),F(3206245775u)},{{F(1036641717u),F(0u)},{F(1036400122u),F(0u)}},{3324,F(1053112776u)},VertexStatus::Valid,0};
 a.vertices[10494]={{F(3212351861u),F(1029544586u),F(3198220100u)},{F(3209537637u),F(3192049056u),F(3205527072u)},{F(1049132724u),F(1060864550u),F(3206569727u)},{{F(1036942486u),F(0u)},{F(1036303796u),F(0u)}},{3308,F(1044115224u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7529];p.patch_id=10284;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130435864u);p.vertex_indices[0]=10364;p.vertex_indices[1]=10365;p.vertex_indices[2]=10493;p.vertex_indices[3]=10494;a.used_patches.push_back(7529);}
 a.vertices[10403]={{F(3212324233u),F(1032954624u),F(1050736452u)},{F(3209629709u),F(3189731826u),F(1058089363u)},{F(3196870136u),F(3208696802u),F(3206090674u)},{{F(1036232646u),F(2147483648u)},{F(3184279266u),F(2147483648u)}},{3296,F(1048844928u)},VertexStatus::Valid,0};
 a.vertices[10404]={{F(3212235016u),F(1032694847u),F(1051291666u)},{F(3209293391u),F(3190062923u),F(1058529861u)},{F(3197032362u),F(3208919658u),F(3205764623u)},{{F(1035924910u),F(2147483648u)},{F(3184359892u),F(2147483648u)}},{3312,F(1051511494u)},VertexStatus::Valid,0};
 a.vertices[10532]={{F(3212351861u),F(1029544586u),F(1050736452u)},{F(3209537637u),F(3192049056u),F(1058043424u)},{F(3196616372u),F(3208348198u),F(3206569727u)},{{F(1036942486u),F(2147483648u)},{F(3183787444u),F(2147483648u)}},{3308,F(1044115224u)},VertexStatus::Valid,0};
 a.vertices[10533]={{F(3212262181u),F(1029022329u),F(1051291666u)},{F(3209199734u),F(3192368336u),F(1058481079u)},{F(3196790039u),F(3208586770u),F(3206245775u)},{{F(1036641717u),F(2147483648u)},{F(3183883770u),F(2147483648u)}},{3324,F(1053112776u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7568];p.patch_id=10323;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130435864u);p.vertex_indices[0]=10403;p.vertex_indices[1]=10404;p.vertex_indices[2]=10532;p.vertex_indices[3]=10533;a.used_patches.push_back(7568);}
 a.vertices[10491]={{F(3212067127u),F(1027886429u),F(3199885741u)},{F(3208471400u),F(3193052469u),F(3206810116u)},{F(1049607771u),F(1061517140u),F(3205641715u)},{{F(1036099754u),F(0u)},{F(1036476167u),F(0u)}},{3359,F(1059373444u)},VertexStatus::Valid,0};
 a.vertices[10492]={{F(3212167305u),F(1028469822u),F(3199330527u)},{F(3208844332u),F(3192702868u),F(3206392610u)},{F(1049464193u),F(1061319901u),F(3205936550u)},{{F(1036360890u),F(0u)},{F(1036456930u),F(0u)}},{3341,F(1057498772u)},VertexStatus::Valid,0};
 a.vertices[10620]={{F(3212087702u),F(1023148776u),F(3199885741u)},{F(3208350656u),F(3195340259u),F(3206750309u)},{F(1049385687u),F(1061212055u),F(3206092242u)},{{F(1036805139u),F(0u)},{F(1036020321u),F(0u)}},{3373,F(1047243308u)},VertexStatus::Valid,0};
 a.vertices[10621]={{F(3212188461u),F(1023866251u),F(3199330527u)},{F(3208725415u),F(3195004244u),F(3206336073u)},{F(1049232309u),F(1061001354u),F(3206385971u)},{{F(1037059634u),F(0u)},{F(1035988812u),F(0u)}},{3354,F(1064809594u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7645];p.patch_id=10410;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130507070u);p.vertex_indices[0]=10491;p.vertex_indices[1]=10492;p.vertex_indices[2]=10620;p.vertex_indices[3]=10621;a.used_patches.push_back(7645);}
 a.vertices[10622]={{F(3212283880u),F(1024421922u),F(3198775314u)},{F(3209082541u),F(3194682550u),F(3205911375u)},{F(1049064500u),F(1060770828u),F(3206692632u)},{{F(1037332535u),F(0u)},{F(1035919211u),F(0u)}},{3337,F(1059984536u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7646];p.patch_id=10411;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130499776u);p.vertex_indices[0]=10492;p.vertex_indices[1]=10493;p.vertex_indices[2]=10621;p.vertex_indices[3]=10622;a.used_patches.push_back(7646);}
 a.vertices[10495]={{F(3212436439u),F(1030037128u),F(3197664886u)},{F(3209858070u),F(3191745209u),F(3205080235u)},{F(1048941374u),F(1060601685u),F(3206908601u)},{{F(1037263316u),F(0u)},{F(1036165663u),F(0u)}},{3292,F(1063819596u)},VertexStatus::Valid,0};
 a.vertices[10623]={{F(3212374064u),F(1024947112u),F(3198220100u)},{F(3209422076u),F(3194375366u),F(3205476825u)},{F(1048880725u),F(1060518370u),F(3207012332u)},{{F(1037623806u),F(0u)},{F(1035809791u),F(0u)}},{3321,F(1051815440u)},VertexStatus::Valid,0};
 a.vertices[10624]={{F(3212459110u),F(1025442379u),F(3197664886u)},{F(3209744043u),F(3194082914u),F(3205033022u)},{F(1048679301u),F(1060241665u),F(3207344999u)},{{F(1037933265u),F(0u)},{F(1035658498u),F(0u)}},{3305,F(1064498082u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7648];p.patch_id=10413;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130486668u);p.vertex_indices[0]=10494;p.vertex_indices[1]=10495;p.vertex_indices[2]=10623;p.vertex_indices[3]=10624;a.used_patches.push_back(7648);}
 a.vertices[10531]={{F(3212436439u),F(1030037128u),F(1050181238u)},{F(3209858070u),F(3191745209u),F(1057596587u)},{F(3196425022u),F(3208085333u),F(3206908601u)},{{F(1037263316u),F(2147483648u)},{F(3183649311u),F(2147483648u)}},{3292,F(1063819596u)},VertexStatus::Valid,0};
 a.vertices[10660]={{F(3212459110u),F(1025442379u),F(1050181238u)},{F(3209744043u),F(3194082914u),F(1057549374u)},{F(3196162949u),F(3207725313u),F(3207344999u)},{{F(1037933265u),F(2147483648u)},{F(3183142146u),F(2147483648u)}},{3305,F(1064498082u)},VertexStatus::Valid,0};
 a.vertices[10661]={{F(3212374064u),F(1024947112u),F(1050736452u)},{F(3209422076u),F(3194375366u),F(1057993177u)},{F(3196364373u),F(3208002018u),F(3207012332u)},{{F(1037623806u),F(2147483648u)},{F(3183293439u),F(2147483648u)}},{3321,F(1051815440u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7685];p.patch_id=10450;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130486667u);p.vertex_indices[0]=10531;p.vertex_indices[1]=10532;p.vertex_indices[2]=10660;p.vertex_indices[3]=10661;a.used_patches.push_back(7685);}
 a.vertices[10534]={{F(3212167305u),F(1028469822u),F(1051846879u)},{F(3208844332u),F(3192702868u),F(1058908962u)},{F(3196947841u),F(3208803549u),F(3205936550u)},{{F(1036360890u),F(2147483648u)},{F(3183940578u),F(2147483648u)}},{3341,F(1057498772u)},VertexStatus::Valid,0};
 a.vertices[10662]={{F(3212283880u),F(1024421922u),F(1051291666u)},{F(3209082541u),F(3194682550u),F(1058427727u)},{F(3196548148u),F(3208254476u),F(3206692632u)},{{F(1037332535u),F(2147483648u)},{F(3183402859u),F(2147483648u)}},{3337,F(1059984536u)},VertexStatus::Valid,0};
 a.vertices[10663]={{F(3212188461u),F(1023866251u),F(1051846879u)},{F(3208725415u),F(3195004244u),F(1058852425u)},{F(3196715957u),F(3208485002u),F(3206385971u)},{{F(1037059634u),F(2147483648u)},{F(3183472460u),F(2147483648u)}},{3354,F(1064809594u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7687];p.patch_id=10452;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130499776u);p.vertex_indices[0]=10533;p.vertex_indices[1]=10534;p.vertex_indices[2]=10662;p.vertex_indices[3]=10663;a.used_patches.push_back(7687);}
 a.vertices[10535]={{F(3212067127u),F(1027886429u),F(1052402093u)},{F(3208471400u),F(3193052469u),F(1059326468u)},{F(3197091419u),F(3209000788u),F(3205641715u)},{{F(1036099754u),F(2147483648u)},{F(3183959815u),F(2147483648u)}},{3359,F(1059373444u)},VertexStatus::Valid,0};
 a.vertices[10664]={{F(3212087702u),F(1023148776u),F(1052402093u)},{F(3208350656u),F(3195340259u),F(1059266661u)},{F(3196869335u),F(3208695703u),F(3206092242u)},{{F(1036805139u),F(2147483648u)},{F(3183503969u),F(2147483648u)}},{3373,F(1047243308u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7688];p.patch_id=10453;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130507070u);p.vertex_indices[0]=10534;p.vertex_indices[1]=10535;p.vertex_indices[2]=10663;p.vertex_indices[3]=10664;a.used_patches.push_back(7688);}
 a.vertices[10750]={{F(3212203979u),F(1015049562u),F(3199330527u)},{F(3208582266u),F(3196685621u),F(3206274508u)},{F(1049001147u),F(1060683796u),F(3206804679u)},{{F(1037737721u),F(0u)},{F(1035515801u),F(0u)}},{3369,F(1060450933u)},VertexStatus::Valid,0};
 a.vertices[10751]={{F(3212299987u),F(1016167776u),F(3198775314u)},{F(3208941169u),F(3196531362u),F(3205853272u)},{F(1048824483u),F(1060441108u),F(3207107016u)},{{F(1038000554u),F(0u)},{F(1035435746u),F(0u)}},{3352,F(1049197174u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7763];p.patch_id=10539;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130558484u);p.vertex_indices[0]=10621;p.vertex_indices[1]=10622;p.vertex_indices[2]=10750;p.vertex_indices[3]=10751;a.used_patches.push_back(7763);}
 a.vertices[10752]={{F(3212390724u),F(1017224591u),F(3198220100u)},{F(3209282381u),F(3196383981u),F(3205422100u)},{F(1048631925u),F(1060176582u),F(3207420720u)},{{F(1038280238u),F(0u)},{F(1035315836u),F(0u)}},{3335,F(1061874774u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7764];p.patch_id=10540;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130551580u);p.vertex_indices[0]=10622;p.vertex_indices[1]=10623;p.vertex_indices[2]=10751;p.vertex_indices[3]=10752;a.used_patches.push_back(7764);}
 a.vertices[10790]={{F(3212390724u),F(1017224591u),F(1050736452u)},{F(3209282381u),F(3196383981u),F(1057938452u)},{F(3196115573u),F(3207660230u),F(3207420720u)},{{F(1038280238u),F(2147483648u)},{F(3182799484u),F(2147483648u)}},{3335,F(1061874774u)},VertexStatus::Valid,0};
 a.vertices[10791]={{F(3212299987u),F(1016167776u),F(1051291666u)},{F(3208941169u),F(3196531362u),F(1058369624u)},{F(3196308131u),F(3207924756u),F(3207107016u)},{{F(1038000554u),F(2147483648u)},{F(3182919394u),F(2147483648u)}},{3352,F(1049197174u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7803];p.patch_id=10579;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130551580u);p.vertex_indices[0]=10661;p.vertex_indices[1]=10662;p.vertex_indices[2]=10790;p.vertex_indices[3]=10791;a.used_patches.push_back(7803);}
 a.vertices[10792]={{F(3212203979u),F(1015049562u),F(1051846879u)},{F(3208582266u),F(3196685621u),F(1058790860u)},{F(3196484795u),F(3208167444u),F(3206804679u)},{{F(1037737721u),F(2147483648u)},{F(3182999449u),F(2147483648u)}},{3369,F(1060450933u)},VertexStatus::Valid,0};
 {auto& p=a.patches[7804];p.patch_id=10580;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3130558484u);p.vertex_indices[0]=10662;p.vertex_indices[1]=10663;p.vertex_indices[2]=10791;p.vertex_indices[3]=10792;a.used_patches.push_back(7804);}
 a.vertices[52399]={{F(3211136324u),F(1033089762u),F(1056310966u)},{F(3197031359u),F(1062059389u),F(3204878341u)},{F(3202164831u),F(3206023241u),F(3207447549u)},{{F(3186581886u),F(1035746630u)},{F(3176819027u),F(1037999975u)}},{13271,F(1052091376u)},VertexStatus::Valid,0};
 a.vertices[52400]={{F(3211102645u),F(1029048238u),F(1056564469u)},{F(3200212551u),F(1062131401u),F(3202955504u)},{F(3202044429u),F(3205890462u),F(3207600076u)},{{F(3190817830u),F(1039535246u)},{F(3180095915u),F(1039557712u)}},{13100,F(1062122797u)},VertexStatus::Valid,0};
 a.vertices[52528]={{F(3211012625u),F(1032263630u),F(1056807505u)},{F(3200400332u),F(1061951600u),F(3203426622u)},{F(3202345675u),F(3206154932u),F(3207272317u)},{{F(3191435592u),F(1039064588u)},{F(3181546863u),F(1039650461u)}},{13104,F(1054592524u)},VertexStatus::Valid,0};
 a.vertices[52529]={{F(3211079669u),F(1032615433u),F(1056546055u)},{F(3201345462u),F(1062019043u),F(3202369628u)},{F(3202288861u),F(3205994710u),F(3207432833u)},{{F(3193180922u),F(1038207442u)},{F(3182320179u),F(1037149010u)}},{13045,F(1064861442u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33340];p.patch_id=51609;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3132790640u);p.vertex_indices[0]=52399;p.vertex_indices[1]=52400;p.vertex_indices[2]=52528;p.vertex_indices[3]=52529;a.used_patches.push_back(33340);}
 a.vertices[52476]={{F(3211102645u),F(1029048238u),F(3204048117u)},{F(3200212551u),F(1062131401u),F(1055471856u)},{F(1054560781u),F(1058406814u),F(3207600076u)},{{F(3190817830u),F(1039535246u)},{F(1032612267u),F(3187041360u)}},{13100,F(1062122797u)},VertexStatus::Valid,0};
 a.vertices[52477]={{F(3211136324u),F(1033089762u),F(3203794614u)},{F(3197031359u),F(1062059389u),F(1057394693u)},{F(1054681183u),F(1058539593u),F(3207447549u)},{{F(3186581886u),F(1035746630u)},{F(1029335379u),F(3185483623u)}},{13271,F(1052091376u)},VertexStatus::Valid,0};
 a.vertices[52605]={{F(3211079669u),F(1032615433u),F(3204029703u)},{F(3201345462u),F(1062019043u),F(1054885980u)},{F(1054805213u),F(1058511062u),F(3207432833u)},{{F(3193180922u),F(1038207442u)},{F(1034836531u),F(3184632658u)}},{13045,F(1064861442u)},VertexStatus::Valid,0};
 a.vertices[52606]={{F(3211012625u),F(1032263630u),F(3204291153u)},{F(3200400332u),F(1061951600u),F(1055942974u)},{F(1054862027u),F(1058671284u),F(3207272317u)},{{F(3191435592u),F(1039064588u)},{F(1034063215u),F(3187134109u)}},{13104,F(1054592524u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33417];p.patch_id=51686;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3132790640u);p.vertex_indices[0]=52476;p.vertex_indices[1]=52477;p.vertex_indices[2]=52605;p.vertex_indices[3]=52606;a.used_patches.push_back(33417);}
 a.vertices[52657]={{F(3211001167u),F(1034800296u),F(1056710158u)},{F(3201331114u),F(1061855504u),F(3202965534u)},{F(3202590308u),F(3206254345u),F(3207100345u)},{{F(3193485451u),F(1036822919u)},{F(3183677477u),F(1036389444u)}},{13060,F(1017062688u)},VertexStatus::Valid,0};
 a.vertices[52658]={{F(3211094839u),F(1036122028u),F(1056281657u)},{F(3201660442u),F(1061955100u),F(3202306456u)},{F(3202557517u),F(3206081605u),F(3207268285u)},{{F(3194644133u),F(1033757915u)},{F(3183982804u),F(1032513060u)}},{13033,F(1034902440u)},VertexStatus::Valid,0};
 a.vertices[52786]={{F(3211020651u),F(1038408298u),F(1056392875u)},{F(3201576010u),F(1061794037u),F(3202957824u)},{F(3202851923u),F(3206338837u),F(3206931558u)},{{F(3194806052u),F(1030735858u)},{F(3185311269u),F(1029479991u)}},{13051,F(1032346664u)},VertexStatus::Valid,0};
 a.vertices[52787]={{F(3211119286u),F(1040209305u),F(1055872767u)},{F(3201588580u),F(1061917723u),F(3202510035u)},{F(3202831677u),F(3206157192u),F(3207107424u)},{{F(3193043960u),F(0u)},{F(3180874247u),F(0u)}},{13040,F(1049406970u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33500];p.patch_id=51865;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=52657;p.vertex_indices[1]=52658;p.vertex_indices[2]=52786;p.vertex_indices[3]=52787;a.used_patches.push_back(33500);}
 a.vertices[52734]={{F(3211094839u),F(1036122028u),F(3203765305u)},{F(3201660442u),F(1061955100u),F(1054822808u)},{F(1055073869u),F(1058597957u),F(3207268285u)},{{F(3194644133u),F(1033757915u)},{F(1036499156u),F(3179996708u)}},{13033,F(1034902440u)},VertexStatus::Valid,0};
 a.vertices[52735]={{F(3211001167u),F(1034800296u),F(3204193806u)},{F(3201331114u),F(1061855504u),F(1055481886u)},{F(1055106660u),F(1058770697u),F(3207100345u)},{{F(3193485451u),F(1036822919u)},{F(1036193829u),F(3183873092u)}},{13060,F(1017062688u)},VertexStatus::Valid,0};
 a.vertices[52863]={{F(3211119286u),F(1040209305u),F(3203356415u)},{F(3201588580u),F(1061917723u),F(1055026387u)},{F(1055348029u),F(1058673544u),F(3207107424u)},{{F(3193043960u),F(0u)},{F(1033390599u),F(0u)}},{13040,F(1049406970u)},VertexStatus::Valid,0};
 a.vertices[52864]={{F(3211020651u),F(1038408298u),F(3203876523u)},{F(3201576010u),F(1061794037u),F(1055474176u)},{F(1055368275u),F(1058855189u),F(3206931558u)},{{F(3194806052u),F(1030735858u)},{F(1037827621u),F(3176963639u)}},{13051,F(1032346664u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33577];p.patch_id=51942;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=52734;p.vertex_indices[1]=52735;p.vertex_indices[2]=52863;p.vertex_indices[3]=52864;a.used_patches.push_back(33577);}
 a.vertices[52915]={{F(3211047503u),F(1041360561u),F(1055936265u)},{F(3201478030u),F(1061750914u),F(3203193723u)},{F(3203115536u),F(3206413499u),F(3206766637u)},{{F(3191498826u),F(0u)},{F(3178942091u),F(0u)}},{13060,F(1023423696u)},VertexStatus::Valid,0};
 a.vertices[52916]={{F(3211140818u),F(1042428028u),F(1055353686u)},{F(3201288324u),F(1061891192u),F(3202879337u)},{F(3203102290u),F(3206224738u),F(3206950702u)},{{F(3190017908u),F(0u)},{F(3172587362u),F(0u)}},{13059,F(1050321718u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33583];p.patch_id=51993;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(955364306u);p.vertex_indices[0]=52786;p.vertex_indices[1]=52787;p.vertex_indices[2]=52915;p.vertex_indices[3]=52916;a.used_patches.push_back(33583);}
 a.vertices[52992]={{F(3211140818u),F(1042428028u),F(3202837334u)},{F(3201288324u),F(1061891192u),F(1055395689u)},{F(1055618642u),F(1058741090u),F(3206950702u)},{{F(3190017908u),F(0u)},{F(1025103714u),F(0u)}},{13059,F(1050321718u)},VertexStatus::Valid,0};
 a.vertices[52993]={{F(3211047503u),F(1041360561u),F(3203419913u)},{F(3201478030u),F(1061750914u),F(1055710075u)},{F(1055631888u),F(1058929851u),F(3206766637u)},{{F(3191498826u),F(0u)},{F(1031458443u),F(0u)}},{13060,F(1023423696u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33660];p.patch_id=52070;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(955364305u);p.vertex_indices[0]=52863;p.vertex_indices[1]=52864;p.vertex_indices[2]=52992;p.vertex_indices[3]=52993;a.used_patches.push_back(33660);}
 a.vertices[53039]={{F(3210705591u),F(1040280499u),F(1057160840u)},{F(3198251747u),F(1060926612u),F(3206107014u)},{F(3203267658u),F(3207300289u),F(3205770948u)},{{F(3189994497u),F(1035882326u)},{F(3186215536u),F(1040963248u)}},{13287,F(1028338528u)},VertexStatus::Valid,0};
 a.vertices[53040]={{F(3210715533u),F(1039906585u),F(1057162319u)},{F(3200419846u),F(1061142523u),F(3205165937u)},{F(3203285749u),F(3207164410u),F(3205918880u)},{{F(3192367918u),F(1035304439u)},{F(3187453526u),F(1038362796u)}},{13163,F(1038390232u)},VertexStatus::Valid,0};
 a.vertices[53168]={{F(3210655557u),F(1040940720u),F(1057191527u)},{F(3200325396u),F(1060903825u),F(3205518323u)},{F(3203445643u),F(3207406979u),F(3205572982u)},{{F(3192326599u),F(1034248976u)},{F(3188333728u),F(1037877607u)}},{13184,F(1059272898u)},VertexStatus::Valid,0};
 a.vertices[53169]={{F(3210726868u),F(1041493861u),F(1057020625u)},{F(3201122364u),F(1061065356u),F(3205018245u)},{F(3203495063u),F(3207256505u),F(3205729470u)},{{F(3193813160u),F(1028400836u)},{F(3188583826u),F(1030907699u)}},{13131,F(1055307072u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33750];p.patch_id=52244;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129773702u);p.vertex_indices[0]=53039;p.vertex_indices[1]=53040;p.vertex_indices[2]=53168;p.vertex_indices[3]=53169;a.used_patches.push_back(33750);}
 a.vertices[53045]={{F(3211153100u),F(1044742446u),F(1054738068u)},{F(3200837922u),F(1061864473u),F(3203361776u)},{F(3203363809u),F(3206286808u),F(3206798102u)},{{F(3188421683u),F(0u)},{F(3165292332u),F(0u)}},{13086,F(1026395792u)},VertexStatus::Valid,0};
 a.vertices[53046]={{F(3211223081u),F(1045946658u),F(1054093372u)},{F(3200415272u),F(1062018250u),F(3203185812u)},{F(3203335448u),F(3206084467u),F(3206997407u)},{{F(3187754281u),F(0u)},{F(3158485285u),F(0u)}},{13096,F(1055241492u)},VertexStatus::Valid,0};
 a.vertices[53174]={{F(3211152333u),F(1047112234u),F(1054031025u)},{F(3200283341u),F(1061829566u),F(3203925327u)},{F(3203612436u),F(3206345743u),F(3206649226u)},{{F(3186859569u),F(0u)},{F(3156993587u),F(0u)}},{13118,F(1032233504u)},VertexStatus::Valid,0};
 a.vertices[53175]={{F(3211204488u),F(1048381303u),F(1053352328u)},{F(3199767299u),F(1061987683u),F(3203798271u)},{F(3203584165u),F(3206139188u),F(3206856452u)},{{F(3185847555u),F(0u)},{F(3146755333u),F(0u)}},{13132,F(1056576633u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33756];p.patch_id=52250;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(970370168u);p.vertex_indices[0]=53045;p.vertex_indices[1]=53046;p.vertex_indices[2]=53174;p.vertex_indices[3]=53175;a.used_patches.push_back(33756);}
 a.vertices[53120]={{F(3211223081u),F(1045946658u),F(3201577020u)},{F(3200415272u),F(1062018250u),F(1055702164u)},{F(1055851800u),F(1058600819u),F(3206997407u)},{{F(3187754281u),F(0u)},{F(1011001637u),F(0u)}},{13096,F(1055241492u)},VertexStatus::Valid,0};
 a.vertices[53121]={{F(3211153100u),F(1044742446u),F(3202221716u)},{F(3200837922u),F(1061864473u),F(1055878128u)},{F(1055880161u),F(1058803160u),F(3206798102u)},{{F(3188421683u),F(0u)},{F(1017808684u),F(0u)}},{13086,F(1026395792u)},VertexStatus::Valid,0};
 a.vertices[53249]={{F(3211204488u),F(1048381303u),F(3200835976u)},{F(3199767299u),F(1061987683u),F(1056314623u)},{F(1056100517u),F(1058655540u),F(3206856452u)},{{F(3185847555u),F(0u)},{F(999271685u),F(0u)}},{13132,F(1056576633u)},VertexStatus::Valid,0};
 a.vertices[53250]={{F(3211152333u),F(1047112234u),F(3201514673u)},{F(3200283341u),F(1061829566u),F(1056441679u)},{F(1056128788u),F(1058862095u),F(3206649226u)},{{F(3186859569u),F(0u)},{F(1009509939u),F(0u)}},{13118,F(1032233504u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33831];p.patch_id=52325;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(970370169u);p.vertex_indices[0]=53120;p.vertex_indices[1]=53121;p.vertex_indices[2]=53249;p.vertex_indices[3]=53250;a.used_patches.push_back(33831);}
 a.vertices[53126]={{F(3210715533u),F(1039906585u),F(3204645967u)},{F(3200419846u),F(1061142523u),F(1057682289u)},{F(1055802101u),F(1059680762u),F(3205918880u)},{{F(3192367918u),F(1035304439u)},{F(1039969878u),F(3185846444u)}},{13163,F(1038390232u)},VertexStatus::Valid,0};
 a.vertices[53127]={{F(3210705591u),F(1040280499u),F(3204644488u)},{F(3198251747u),F(1060926612u),F(1058623366u)},{F(1055784010u),F(1059816641u),F(3205770948u)},{{F(3189994497u),F(1035882326u)},{F(1038731888u),F(3188446896u)}},{13287,F(1028338528u)},VertexStatus::Valid,0};
 a.vertices[53255]={{F(3210726868u),F(1041493861u),F(3204504273u)},{F(3201122364u),F(1061065356u),F(1057534597u)},{F(1056011415u),F(1059772857u),F(3205729470u)},{{F(3193813160u),F(1028400836u)},{F(1041100178u),F(3178391347u)}},{13131,F(1055307072u)},VertexStatus::Valid,0};
 a.vertices[53256]={{F(3210655557u),F(1040940720u),F(3204675175u)},{F(3200325396u),F(1060903825u),F(1058034675u)},{F(1055961995u),F(1059923331u),F(3205572982u)},{{F(3192326599u),F(1034248976u)},{F(1040850080u),F(3185361255u)}},{13184,F(1059272898u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33837];p.patch_id=52331;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3129773701u);p.vertex_indices[0]=53126;p.vertex_indices[1]=53127;p.vertex_indices[2]=53255;p.vertex_indices[3]=53256;a.used_patches.push_back(33837);}
 a.vertices[53297]={{F(3210669364u),F(1042354564u),F(1057039081u)},{F(3201063931u),F(1060827865u),F(3205361498u)},{F(3203640330u),F(3207499346u),F(3205379183u)},{{F(3193706374u),F(1024244838u)},{F(3189368917u),F(1026661627u)}},{13151,F(1053763802u)},VertexStatus::Valid,0};
 a.vertices[53298]={{F(3210761335u),F(1043276773u),F(1056600554u)},{F(3201329846u),F(1060987273u),F(3205048226u)},{F(3203706369u),F(3207339988u),F(3205543630u)},{{F(3190522493u),F(0u)},{F(3181402308u),F(0u)}},{13126,F(1019799360u)},VertexStatus::Valid,0};
 a.vertices[53426]={{F(3210707348u),F(1044088515u),F(1056612711u)},{F(3201308424u),F(1060743563u),F(3205383320u)},{F(3203834995u),F(3207584199u),F(3205188796u)},{{F(3189823856u),F(0u)},{F(3181055269u),F(0u)}},{13144,F(1055719510u)},VertexStatus::Valid,0};
 a.vertices[53427]={{F(3210802933u),F(1045230890u),F(1055986137u)},{F(3201271626u),F(1060906293u),F(3205180772u)},{F(3203910919u),F(3207418510u),F(3205361145u)},{{F(3188557511u),F(0u)},{F(3175047253u),F(0u)}},{13134,F(1063527014u)},VertexStatus::Valid,0};
 {auto& p=a.patches[33930];p.patch_id=52500;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=53297;p.vertex_indices[1]=53298;p.vertex_indices[2]=53426;p.vertex_indices[3]=53427;a.used_patches.push_back(33930);}
 a.vertices[53384]={{F(3210761335u),F(1043276773u),F(3204084202u)},{F(3201329846u),F(1060987273u),F(1057564578u)},{F(1056222721u),F(1059856340u),F(3205543630u)},{{F(3190522493u),F(0u)},{F(1033918660u),F(0u)}},{13126,F(1019799360u)},VertexStatus::Valid,0};
 a.vertices[53385]={{F(3210669364u),F(1042354564u),F(3204522729u)},{F(3201063931u),F(1060827865u),F(1057877850u)},{F(1056156682u),F(1060015698u),F(3205379183u)},{{F(3193706374u),F(1024244838u)},{F(1041885269u),F(3174145275u)}},{13151,F(1053763802u)},VertexStatus::Valid,0};
 a.vertices[53513]={{F(3210802933u),F(1045230890u),F(3203469785u)},{F(3201271626u),F(1060906293u),F(1057697124u)},{F(1056427271u),F(1059934862u),F(3205361145u)},{{F(3188557511u),F(0u)},{F(1027563605u),F(0u)}},{13134,F(1063527014u)},VertexStatus::Valid,0};
 a.vertices[53514]={{F(3210707348u),F(1044088515u),F(3204096359u)},{F(3201308424u),F(1060743563u),F(1057899672u)},{F(1056351347u),F(1060100551u),F(3205188796u)},{{F(3189823856u),F(0u)},{F(1033571621u),F(0u)}},{13144,F(1055719510u)},VertexStatus::Valid,0};
 {auto& p=a.patches[34017];p.patch_id=52587;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=53384;p.vertex_indices[1]=53385;p.vertex_indices[2]=53513;p.vertex_indices[3]=53514;a.used_patches.push_back(34017);}
 a.vertices[53555]={{F(3210753693u),F(1045983698u),F(1055970891u)},{F(3201289418u),F(1060652414u),F(3205507935u)},{F(3204021547u),F(3207665025u),F(3205001428u)},{{F(3188104213u),F(0u)},{F(3175061375u),F(0u)}},{13151,F(1061817174u)},VertexStatus::Valid,0};
 a.vertices[53556]={{F(3210843513u),F(1047274267u),F(1055256047u)},{F(3201050090u),F(1060817503u),F(3205380101u)},{F(3204103399u),F(3207494712u),F(3205181526u)},{{F(3186869194u),F(0u)},{F(3169718212u),F(0u)}},{13152,F(1059060549u)},VertexStatus::Valid,0};
 {auto& p=a.patches[34022];p.patch_id=52628;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(953681430u);p.vertex_indices[0]=53426;p.vertex_indices[1]=53427;p.vertex_indices[2]=53555;p.vertex_indices[3]=53556;a.used_patches.push_back(34022);}
 a.vertices[53642]={{F(3210843513u),F(1047274267u),F(3202739695u)},{F(3201050090u),F(1060817503u),F(1057896453u)},{F(1056619751u),F(1060011064u),F(3205181526u)},{{F(3186869194u),F(0u)},{F(1022234564u),F(0u)}},{13152,F(1059060549u)},VertexStatus::Valid,0};
 a.vertices[53643]={{F(3210753693u),F(1045983698u),F(3203454539u)},{F(3201289418u),F(1060652414u),F(1058024287u)},{F(1056537899u),F(1060181377u),F(3205001428u)},{{F(3188104213u),F(0u)},{F(1027577727u),F(0u)}},{13151,F(1061817174u)},VertexStatus::Valid,0};
 {auto& p=a.patches[34109];p.patch_id=52715;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(953681430u);p.vertex_indices[0]=53513;p.vertex_indices[1]=53514;p.vertex_indices[2]=53642;p.vertex_indices[3]=53643;a.used_patches.push_back(34109);}
 a.vertices[53557]={{F(3210920440u),F(1048606366u),F(1054501204u)},{F(3200659662u),F(1060981305u),F(3205301256u)},{F(3204171192u),F(3207315753u),F(3205371541u)},{{F(3185563840u),F(0u)},{F(3163340200u),F(0u)}},{13161,F(1010363712u)},VertexStatus::Valid,0};
 a.vertices[53558]={{F(3210980799u),F(1049307493u),F(1053723838u)},{F(3200161437u),F(1061144205u),F(3205252728u)},{F(3204222267u),F(3207128833u),F(3205570818u)},{{F(3184503911u),F(0u)},{F(3155284576u),F(0u)}},{13174,F(1056871290u)},VertexStatus::Valid,0};
 a.vertices[53686]={{F(3210938631u),F(1049687870u),F(1053603241u)},{F(3200216909u),F(1060879623u),F(3205585728u)},{F(3204349589u),F(3207389317u),F(3205201315u)},{{F(3183869374u),F(0u)},{F(3155512840u),F(0u)}},{13190,F(1026907168u)},VertexStatus::Valid,0};
 a.vertices[53687]={{F(3210979093u),F(1050420735u),F(1052771590u)},{F(3199626273u),F(1061040899u),F(3205565930u)},{F(3204401066u),F(3207200368u),F(3205407689u)},{{F(3183015349u),F(0u)},{F(3142340845u),F(0u)}},{13207,F(1060076938u)},VertexStatus::Valid,0};
 {auto& p=a.patches[34117];p.patch_id=52758;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(971065652u);p.vertex_indices[0]=53557;p.vertex_indices[1]=53558;p.vertex_indices[2]=53686;p.vertex_indices[3]=53687;a.used_patches.push_back(34117);}
 a.vertices[53640]={{F(3210980799u),F(1049307493u),F(3201207486u)},{F(3200161437u),F(1061144205u),F(1057769080u)},{F(1056738619u),F(1059645185u),F(3205570818u)},{{F(3184503911u),F(0u)},{F(1007800928u),F(0u)}},{13174,F(1056871290u)},VertexStatus::Valid,0};
 a.vertices[53641]={{F(3210920440u),F(1048606366u),F(3201984852u)},{F(3200659662u),F(1060981305u),F(1057817608u)},{F(1056687544u),F(1059832105u),F(3205371541u)},{{F(3185563840u),F(0u)},{F(1015856552u),F(0u)}},{13161,F(1010363712u)},VertexStatus::Valid,0};
 a.vertices[53769]={{F(3210979093u),F(1050420735u),F(3200255238u)},{F(3199626273u),F(1061040899u),F(1058082282u)},{F(1056917418u),F(1059716720u),F(3205407689u)},{{F(3183015349u),F(0u)},{F(994857197u),F(0u)}},{13207,F(1060076938u)},VertexStatus::Valid,0};
 a.vertices[53770]={{F(3210938631u),F(1049687870u),F(3201086889u)},{F(3200216909u),F(1060879623u),F(1058102080u)},{F(1056865941u),F(1059905669u),F(3205201315u)},{{F(3183869374u),F(0u)},{F(1008029192u),F(0u)}},{13190,F(1026907168u)},VertexStatus::Valid,0};
 {auto& p=a.patches[34200];p.patch_id=52841;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(971065652u);p.vertex_indices[0]=53640;p.vertex_indices[1]=53641;p.vertex_indices[2]=53769;p.vertex_indices[3]=53770;a.used_patches.push_back(34200);}
 a.vertices[54578]={{F(3210435556u),F(1048297368u),F(1056464204u)},{F(3200112591u),F(1057723342u),F(3208673255u)},{F(3203558228u),F(3209605666u),F(3199482845u)},{{F(3185785458u),F(0u)},{F(3189746199u),F(0u)}},{13382,F(1054005604u)},VertexStatus::Valid,0};
 a.vertices[54579]={{F(3210472048u),F(1048696193u),F(1056174685u)},{F(3201676352u),F(1058097832u),F(3207966475u)},{F(3203721987u),F(3209502302u),F(3199731191u)},{{F(3184591357u),F(0u)},{F(3183856090u),F(0u)}},{13288,F(1052111675u)},VertexStatus::Valid,0};
 a.vertices[54707]={{F(3210450753u),F(1048756555u),F(1056213359u)},{F(3201454678u),F(1057695114u),F(3208340510u)},{F(3203563790u),F(3209703665u),F(3199014926u)},{{F(3184272434u),F(0u)},{F(3186217131u),F(0u)}},{13320,F(1052883208u)},VertexStatus::Valid,0};
 a.vertices[54708]={{F(3210543440u),F(1049225226u),F(1055563591u)},{F(3202200108u),F(1057963038u),F(3207916504u)},{F(3203731875u),F(3209600260u),F(3199270906u)},{{F(3182994675u),F(0u)},{F(3180501529u),F(0u)}},{13270,F(1042517559u)},VertexStatus::Valid,0};
 {auto& p=a.patches[34910];p.patch_id=53771;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3128763700u);p.vertex_indices[0]=54578;p.vertex_indices[1]=54579;p.vertex_indices[2]=54707;p.vertex_indices[3]=54708;a.used_patches.push_back(34910);}
 a.vertices[54683]={{F(3210472048u),F(1048696193u),F(3203658333u)},{F(3201676352u),F(1058097832u),F(1060482827u)},{F(1056238339u),F(1062018654u),F(3199731191u)},{{F(3184591357u),F(0u)},{F(1036372442u),F(0u)}},{13288,F(1052111675u)},VertexStatus::Valid,0};
 a.vertices[54684]={{F(3210435556u),F(1048297368u),F(3203947852u)},{F(3200112591u),F(1057723342u),F(1061189607u)},{F(1056074580u),F(1062122018u),F(3199482845u)},{{F(3185785458u),F(0u)},{F(1042262551u),F(0u)}},{13382,F(1054005604u)},VertexStatus::Valid,0};
 a.vertices[54812]={{F(3210543440u),F(1049225226u),F(3203047239u)},{F(3202200108u),F(1057963038u),F(1060432856u)},{F(1056248227u),F(1062116612u),F(3199270906u)},{{F(3182994675u),F(0u)},{F(1033017881u),F(0u)}},{13270,F(1042517559u)},VertexStatus::Valid,0};
 a.vertices[54813]={{F(3210450753u),F(1048756555u),F(3203697007u)},{F(3201454678u),F(1057695114u),F(1060856862u)},{F(1056080142u),F(1062220017u),F(3199014926u)},{{F(3184272434u),F(0u)},{F(1038733483u),F(0u)}},{13320,F(1052883208u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35015];p.patch_id=53876;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3128763700u);p.vertex_indices[0]=54683;p.vertex_indices[1]=54684;p.vertex_indices[2]=54812;p.vertex_indices[3]=54813;a.used_patches.push_back(35015);}
 a.vertices[54835]={{F(3210446759u),F(1048852356u),F(1056166958u)},{F(3200943078u),F(1057241368u),F(3208801298u)},{F(3203387891u),F(3209896118u),F(3198305368u)},{{F(3183767665u),F(0u)},{F(3188409871u),F(0u)}},{13367,F(1060294146u)},VertexStatus::Valid,0};
 a.vertices[54836]={{F(3210509477u),F(1049193029u),F(1055712823u)},{F(3202193163u),F(1057580543u),F(3208211611u)},{F(3203555042u),F(3209801465u),F(3198545239u)},{{F(3182819854u),F(0u)},{F(3182434219u),F(0u)}},{13290,F(1061031348u)},VertexStatus::Valid,0};
 a.vertices[54964]={{F(3210479919u),F(1049169507u),F(1055838228u)},{F(3202020204u),F(1057170019u),F(3208555375u)},{F(3203361739u),F(3209992958u),F(3197827478u)},{{F(3182615298u),F(0u)},{F(3184944960u),F(0u)}},{13320,F(1027758784u)},VertexStatus::Valid,0};
 a.vertices[54965]={{F(3210590922u),F(1049670600u),F(1055062636u)},{F(3202643450u),F(1057422985u),F(3208190706u)},{F(3203529307u),F(3209900123u),F(3198072008u)},{{F(3181509516u),F(0u)},{F(3179720753u),F(0u)}},{13277,F(1032910360u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35124];p.patch_id=54026;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3125957074u);p.vertex_indices[0]=54835;p.vertex_indices[1]=54836;p.vertex_indices[2]=54964;p.vertex_indices[3]=54965;a.used_patches.push_back(35124);}
 a.vertices[54837]={{F(3210630448u),F(1049777775u),F(1054825981u)},{F(3202506525u),F(1057796055u),F(3207952616u)},{F(3203724026u),F(3209699518u),F(3198807348u)},{{F(3181622737u),F(0u)},{F(3176965002u),F(0u)}},{13264,F(1044025320u)},VertexStatus::Valid,0};
 a.vertices[54966]={{F(3210725547u),F(1050324011u),F(1053994951u)},{F(3202680912u),F(1057604372u),F(3208044614u)},{F(3203696480u),F(3209800877u),F(3198338737u)},{{F(3180438703u),F(0u)},{F(3173483672u),F(0u)}},{13265,F(1062517556u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35125];p.patch_id=54027;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3111131688u);p.vertex_indices[0]=54836;p.vertex_indices[1]=54837;p.vertex_indices[2]=54965;p.vertex_indices[3]=54966;a.used_patches.push_back(35125);}
 a.vertices[54838]={{F(3210756977u),F(1050475291u),F(1053729174u)},{F(3202389008u),F(1057965324u),F(3207856297u)},{F(3203890725u),F(3209590807u),F(3199092629u)},{{F(3180568411u),F(0u)},{F(3171431260u),F(0u)}},{13260,F(1062295540u)},VertexStatus::Valid,0};
 a.vertices[54967]={{F(3210850652u),F(1051054434u),F(1052764917u)},{F(3202400239u),F(1057752592u),F(3208019072u)},{F(3203859777u),F(3209695683u),F(3198628179u)},{{F(3179527800u),F(0u)},{F(3167396086u),F(0u)}},{13271,F(1057029825u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35126];p.patch_id=54028;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=54837;p.vertex_indices[1]=54838;p.vertex_indices[2]=54966;p.vertex_indices[3]=54967;a.used_patches.push_back(35126);}
 a.vertices[54940]={{F(3210756977u),F(1050475291u),F(3201212822u)},{F(3202389008u),F(1057965324u),F(1060372649u)},{F(1056407077u),F(1062107159u),F(3199092629u)},{{F(3180568411u),F(0u)},{F(1023947612u),F(0u)}},{13260,F(1062295540u)},VertexStatus::Valid,0};
 a.vertices[54941]={{F(3210630448u),F(1049777775u),F(3202309629u)},{F(3202506525u),F(1057796055u),F(1060468968u)},{F(1056240378u),F(1062215870u),F(3198807348u)},{{F(3181622737u),F(0u)},{F(1029481354u),F(0u)}},{13264,F(1044025320u)},VertexStatus::Valid,0};
 a.vertices[55069]={{F(3210850652u),F(1051054434u),F(3200248565u)},{F(3202400239u),F(1057752592u),F(1060535424u)},{F(1056376129u),F(1062212035u),F(3198628179u)},{{F(3179527800u),F(0u)},{F(1019912438u),F(0u)}},{13271,F(1057029825u)},VertexStatus::Valid,0};
 a.vertices[55070]={{F(3210725547u),F(1050324011u),F(3201478599u)},{F(3202680912u),F(1057604372u),F(1060560966u)},{F(1056212832u),F(1062317229u),F(3198338737u)},{{F(3180438703u),F(0u)},{F(1026000024u),F(0u)}},{13265,F(1062517556u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35229];p.patch_id=54131;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=54940;p.vertex_indices[1]=54941;p.vertex_indices[2]=55069;p.vertex_indices[3]=55070;a.used_patches.push_back(35229);}
 a.vertices[54942]={{F(3210509477u),F(1049193029u),F(3203196471u)},{F(3202193163u),F(1057580543u),F(1060727963u)},{F(1056071394u),F(1062317817u),F(3198545239u)},{{F(3182819854u),F(0u)},{F(1034950571u),F(0u)}},{13290,F(1061031348u)},VertexStatus::Valid,0};
 a.vertices[55071]={{F(3210590922u),F(1049670600u),F(3202546284u)},{F(3202643450u),F(1057422985u),F(1060707058u)},{F(1056045659u),F(1062416475u),F(3198072008u)},{{F(3181509516u),F(0u)},{F(1032237105u),F(0u)}},{13277,F(1032910360u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35230];p.patch_id=54132;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3111131689u);p.vertex_indices[0]=54941;p.vertex_indices[1]=54942;p.vertex_indices[2]=55070;p.vertex_indices[3]=55071;a.used_patches.push_back(35230);}
 a.vertices[54943]={{F(3210446759u),F(1048852356u),F(3203650606u)},{F(3200943078u),F(1057241368u),F(1061317650u)},{F(1055904243u),F(1062412470u),F(3198305368u)},{{F(3183767665u),F(0u)},{F(1040926223u),F(0u)}},{13367,F(1060294146u)},VertexStatus::Valid,0};
 a.vertices[55072]={{F(3210479919u),F(1049169507u),F(3203321876u)},{F(3202020204u),F(1057170019u),F(1061071727u)},{F(1055878091u),F(1062509310u),F(3197827478u)},{{F(3182615298u),F(0u)},{F(1037461312u),F(0u)}},{13320,F(1027758784u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35231];p.patch_id=54133;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3125957074u);p.vertex_indices[0]=54942;p.vertex_indices[1]=54943;p.vertex_indices[2]=55071;p.vertex_indices[3]=55072;a.used_patches.push_back(35231);}
 a.vertices[55094]={{F(3210684737u),F(1050151670u),F(1054309959u)},{F(3202927056u),F(1057236203u),F(3208237060u)},{F(3203484736u),F(3210000323u),F(3197593603u)},{{F(3180353382u),F(0u)},{F(3176156889u),F(0u)}},{13272,F(1063522101u)},VertexStatus::Valid,0};
 a.vertices[55095]={{F(3210824137u),F(1050844480u),F(1053092214u)},{F(3202770410u),F(1057391248u),F(3208174460u)},{F(3203647887u),F(3209904843u),F(3197863316u)},{{F(3179411653u),F(0u)},{F(3170876085u),F(0u)}},{13272,F(1057180533u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35234];p.patch_id=54155;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=54965;p.vertex_indices[1]=54966;p.vertex_indices[2]=55094;p.vertex_indices[3]=55095;a.used_patches.push_back(35234);}
 a.vertices[54968]={{F(3210950134u),F(1051826487u),F(1051432436u)},{F(3201915156u),F(1057881399u),F(3208066347u)},{F(3204016445u),F(3209584962u),F(3198940080u)},{{F(3178258337u),F(0u)},{F(3160768087u),F(0u)}},{13287,F(1057279908u)},VertexStatus::Valid,0};
 a.vertices[54969]={{F(3211015025u),F(1052620190u),F(1050032185u)},{F(3201285409u),F(1057996785u),F(3208158659u)},{F(3204164334u),F(3209469050u),F(3199273748u)},{{F(3176999448u),F(0u)},{F(3151919351u),F(0u)}},{13310,F(1050660206u)},VertexStatus::Valid,0};
 a.vertices[55097]={{F(3211029313u),F(1052363749u),F(1050306220u)},{F(3201788472u),F(1057632232u),F(3208292154u)},{F(3203956738u),F(3209698661u),F(3198468501u)},{{F(3176658063u),F(0u)},{F(3156210217u),F(0u)}},{13305,F(1064919620u)},VertexStatus::Valid,0};
 a.vertices[55098]={{F(3211075816u),F(1053149367u),F(1048811783u)},{F(3201087747u),F(1057731668u),F(3208414677u)},{F(3204098306u),F(3209588589u),F(3198802790u)},{{F(3175570515u),F(0u)},{F(3145927357u),F(0u)}},{13332,F(1052878976u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35237];p.patch_id=54158;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(972146116u);p.vertex_indices[0]=54968;p.vertex_indices[1]=54969;p.vertex_indices[2]=55097;p.vertex_indices[3]=55098;a.used_patches.push_back(35237);}
 a.vertices[55067]={{F(3211015025u),F(1052620190u),F(3197515833u)},{F(3201285409u),F(1057996785u),F(1060675011u)},{F(1056680686u),F(1061985402u),F(3199273748u)},{{F(3176999448u),F(0u)},{F(1004435703u),F(0u)}},{13310,F(1050660206u)},VertexStatus::Valid,0};
 a.vertices[55068]={{F(3210950134u),F(1051826487u),F(3198916084u)},{F(3201915156u),F(1057881399u),F(1060582699u)},{F(1056532797u),F(1062101314u),F(3198940080u)},{{F(3178258337u),F(0u)},{F(1013284439u),F(0u)}},{13287,F(1057279908u)},VertexStatus::Valid,0};
 a.vertices[55196]={{F(3211075816u),F(1053149367u),F(3196295431u)},{F(3201087747u),F(1057731668u),F(1060931029u)},{F(1056614658u),F(1062104941u),F(3198802790u)},{{F(3175570515u),F(0u)},{F(998443709u),F(0u)}},{13332,F(1052878976u)},VertexStatus::Valid,0};
 a.vertices[55197]={{F(3211029313u),F(1052363749u),F(3197789868u)},{F(3201788472u),F(1057632232u),F(1060808506u)},{F(1056473090u),F(1062215013u),F(3198468501u)},{{F(3176658063u),F(0u)},{F(1008726569u),F(0u)}},{13305,F(1064919620u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35336];p.patch_id=54257;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(972146116u);p.vertex_indices[0]=55067;p.vertex_indices[1]=55068;p.vertex_indices[2]=55196;p.vertex_indices[3]=55197;a.used_patches.push_back(35336);}
 a.vertices[55199]={{F(3210824137u),F(1050844480u),F(3200575862u)},{F(3202770410u),F(1057391248u),F(1060690812u)},{F(1056164239u),F(1062421195u),F(3197863316u)},{{F(3179411653u),F(0u)},{F(1023392437u),F(0u)}},{13272,F(1057180533u)},VertexStatus::Valid,0};
 a.vertices[55200]={{F(3210684737u),F(1050151670u),F(3201793607u)},{F(3202927056u),F(1057236203u),F(1060753412u)},{F(1056001088u),F(1062516675u),F(3197593603u)},{{F(3180353382u),F(0u)},{F(1028673241u),F(0u)}},{13272,F(1063522101u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35339];p.patch_id=54260;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=55070;p.vertex_indices[1]=55071;p.vertex_indices[2]=55199;p.vertex_indices[3]=55200;a.used_patches.push_back(35339);}
 a.vertices[55096]={{F(3210943062u),F(1051590336u),F(1051742544u)},{F(3202364352u),F(1057520199u),F(3208204728u)},{F(3203805993u),F(3209804173u),F(3198155088u)},{{F(3177965065u),F(0u)},{F(3163963907u),F(0u)}},{13285,F(1038203560u)},VertexStatus::Valid,0};
 a.vertices[55224]={{F(3210923252u),F(1051325374u),F(1052134332u)},{F(3202804172u),F(1057158347u),F(3208330158u)},{F(3203577337u),F(3210011686u),F(3197379391u)},{{F(3177748425u),F(0u)},{F(3166945629u),F(0u)}},{13282,F(1058596340u)},VertexStatus::Valid,0};
 a.vertices[55225]={{F(3211032595u),F(1052072717u),F(1050674631u)},{F(3202298074u),F(1057268666u),F(3208405700u)},{F(3203728677u),F(3209916432u),F(3197671560u)},{{F(3176382052u),F(0u)},{F(3160235491u),F(0u)}},{13300,F(1059973125u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35345];p.patch_id=54284;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(965713689u);p.vertex_indices[0]=55095;p.vertex_indices[1]=55096;p.vertex_indices[2]=55224;p.vertex_indices[3]=55225;a.used_patches.push_back(35345);}
 a.vertices[55198]={{F(3210943062u),F(1051590336u),F(3199226192u)},{F(3202364352u),F(1057520199u),F(1060721080u)},{F(1056322345u),F(1062320525u),F(3198155088u)},{{F(3177965065u),F(0u)},{F(1016480259u),F(0u)}},{13285,F(1038203560u)},VertexStatus::Valid,0};
 a.vertices[55327]={{F(3211032595u),F(1052072717u),F(3198158279u)},{F(3202298074u),F(1057268666u),F(1060922052u)},{F(1056245029u),F(1062432784u),F(3197671560u)},{{F(3176382052u),F(0u)},{F(1012751843u),F(0u)}},{13300,F(1059973125u)},VertexStatus::Valid,0};
 a.vertices[55328]={{F(3210923252u),F(1051325374u),F(3199617980u)},{F(3202804172u),F(1057158347u),F(1060846510u)},{F(1056093689u),F(1062528038u),F(3197379391u)},{{F(3177748425u),F(0u)},{F(1019461981u),F(0u)}},{13282,F(1058596340u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35448];p.patch_id=54387;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(965713688u);p.vertex_indices[0]=55198;p.vertex_indices[1]=55199;p.vertex_indices[2]=55327;p.vertex_indices[3]=55328;a.used_patches.push_back(35448);}
 a.vertices[55226]={{F(3211103391u),F(1052836448u),F(1049140216u)},{F(3201646721u),F(1057363554u),F(3208526563u)},{F(3203872135u),F(3209816995u),F(3197984432u)},{{F(3175249562u),F(0u)},{F(3151488168u),F(0u)}},{13325,F(1057189432u)},VertexStatus::Valid,0};
 a.vertices[55354]={{F(3211118232u),F(1052493485u),F(1049573247u)},{F(3202213663u),F(1056998509u),F(3208616096u)},{F(3203627416u),F(3210032402u),F(3197175916u)},{{F(3174986120u),F(0u)},{F(3156305120u),F(0u)}},{13317,F(1056909410u)},VertexStatus::Valid,0};
 a.vertices[55355]={{F(3211171953u),F(1053237618u),F(1047315799u)},{F(3201498432u),F(1057075758u),F(3208765156u)},{F(3203762310u),F(3209939803u),F(3197486134u)},{{F(3174000968u),F(0u)},{F(3146416309u),F(0u)}},{13345,F(1059278360u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35457];p.patch_id=54413;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(972703260u);p.vertex_indices[0]=55225;p.vertex_indices[1]=55226;p.vertex_indices[2]=55354;p.vertex_indices[3]=55355;a.used_patches.push_back(35457);}
 a.vertices[55326]={{F(3211103391u),F(1052836448u),F(3196623864u)},{F(3201646721u),F(1057363554u),F(1061042915u)},{F(1056388487u),F(1062333347u),F(3197984432u)},{{F(3175249562u),F(0u)},{F(1004004520u),F(0u)}},{13325,F(1057189432u)},VertexStatus::Valid,0};
 a.vertices[55455]={{F(3211171953u),F(1053237618u),F(3194799447u)},{F(3201498432u),F(1057075758u),F(1061281508u)},{F(1056278662u),F(1062456155u),F(3197486134u)},{{F(3174000968u),F(0u)},{F(998932661u),F(0u)}},{13345,F(1059278360u)},VertexStatus::Valid,0};
 a.vertices[55456]={{F(3211118232u),F(1052493485u),F(3197056895u)},{F(3202213663u),F(1056998509u),F(1061132448u)},{F(1056143768u),F(1062548754u),F(3197175916u)},{{F(3174986120u),F(0u)},{F(1008821472u),F(0u)}},{13317,F(1056909410u)},VertexStatus::Valid,0};
 {auto& p=a.patches[35558];p.patch_id=54514;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(972703259u);p.vertex_indices[0]=55326;p.vertex_indices[1]=55327;p.vertex_indices[2]=55455;p.vertex_indices[3]=55456;a.used_patches.push_back(35558);}
 a.vertices[60902]={{F(3211383121u),F(1051372472u),F(1049460680u)},{F(3208631666u),F(3194066626u),F(3206536220u)},{F(3192061074u),F(3210010095u),F(1057331235u)},{{F(1031473088u),F(0u)},{F(3164563997u),F(0u)}},{13189,F(1036492340u)},VertexStatus::Valid,0};
 a.vertices[60903]={{F(3211703109u),F(1050692857u),F(1047714331u)},{F(3208472665u),F(3195310736u),F(3206611206u)},{F(3191708842u),F(3209925184u),F(1057496253u)},{{F(1030517649u),F(0u)},{F(3160174691u),F(0u)}},{13211,F(1061929493u)},VertexStatus::Valid,0};
 a.vertices[61031]={{F(3211137989u),F(1051912224u),F(1050269126u)},{F(3208798075u),F(3194280116u),F(3206313881u)},{F(3191870026u),F(3209887052u),F(1057541205u)},{{F(1032374609u),F(0u)},{F(3167021705u),F(0u)}},{13172,F(1051882804u)},VertexStatus::Valid,0};
 a.vertices[61032]={{F(3211472222u),F(1051266790u),F(1049004478u)},{F(3208662374u),F(3195463466u),F(3206368971u)},{F(3191523966u),F(3209797248u),F(1057707669u)},{{F(1031829273u),F(0u)},{F(3163167996u),F(0u)}},{13192,F(1060840294u)},VertexStatus::Valid,0};
 {auto& p=a.patches[40688];p.patch_id=60046;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(968264599u);p.vertex_indices[0]=60902;p.vertex_indices[1]=60903;p.vertex_indices[2]=61031;p.vertex_indices[3]=61032;a.used_patches.push_back(40688);}
 a.vertices[61001]={{F(3211703109u),F(1050692857u),F(3195197979u)},{F(3208472665u),F(3195310736u),F(1059127558u)},{F(1044225194u),F(1062441536u),F(1057496253u)},{{F(1030517649u),F(0u)},{F(1012691043u),F(0u)}},{13211,F(1061929493u)},VertexStatus::Valid,0};
 a.vertices[61002]={{F(3211383121u),F(1051372472u),F(3196944328u)},{F(3208631666u),F(3194066626u),F(1059052572u)},{F(1044577426u),F(1062526447u),F(1057331235u)},{{F(1031473088u),F(0u)},{F(1017080349u),F(0u)}},{13189,F(1036492340u)},VertexStatus::Valid,0};
 a.vertices[61130]={{F(3211472222u),F(1051266790u),F(3196488126u)},{F(3208662374u),F(3195463466u),F(1058885323u)},{F(1044040318u),F(1062313600u),F(1057707669u)},{{F(1031829273u),F(0u)},{F(1015684348u),F(0u)}},{13192,F(1060840294u)},VertexStatus::Valid,0};
 a.vertices[61131]={{F(3211137989u),F(1051912224u),F(3197752774u)},{F(3208798075u),F(3194280116u),F(1058830233u)},{F(1044386378u),F(1062403404u),F(1057541205u)},{{F(1032374609u),F(0u)},{F(1019538057u),F(0u)}},{13172,F(1051882804u)},VertexStatus::Valid,0};
 {auto& p=a.patches[40787];p.patch_id=60145;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(968264598u);p.vertex_indices[0]=61001;p.vertex_indices[1]=61002;p.vertex_indices[2]=61130;p.vertex_indices[3]=61131;a.used_patches.push_back(40787);}
 a.vertices[61029]={{F(3210400015u),F(1053014123u),F(1052744866u)},{F(3208959925u),F(3192287929u),F(3206281508u)},{F(3192546802u),F(3210055503u),F(1057212811u)},{{F(1033789577u),F(0u)},{F(3175620116u),F(0u)}},{13143,F(1062677489u)},VertexStatus::Valid,0};
 a.vertices[61030]={{F(3210778655u),F(1052496443u),F(1051521344u)},{F(3208900786u),F(3193208788u),F(3206280152u)},{F(3192210516u),F(3209973075u),F(1057376315u)},{{F(1033019107u),F(0u)},{F(3171681836u),F(0u)}},{13155,F(1056935908u)},VertexStatus::Valid,0};
 a.vertices[61158]={{F(3210130998u),F(1053480497u),F(1053396127u)},{F(3209044897u),F(3192713480u),F(3206138523u)},{F(3192351702u),F(3209943002u),F(1057410709u)},{{F(1034826600u),F(0u)},{F(3178554871u),F(0u)}},{13136,F(1058083680u)},VertexStatus::Valid,0};
 a.vertices[61159]={{F(3210510398u),F(1053003800u),F(1052240521u)},{F(3209021109u),F(3193516742u),F(3206100169u)},{F(3192019465u),F(3209855711u),F(1057576229u)},{{F(1033937054u),F(0u)},{F(3173759671u),F(0u)}},{13143,F(1060408307u)},VertexStatus::Valid,0};
 {auto& p=a.patches[40798];p.patch_id=60172;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(958374564u);p.vertex_indices[0]=61029;p.vertex_indices[1]=61030;p.vertex_indices[2]=61158;p.vertex_indices[3]=61159;a.used_patches.push_back(40798);}
 a.vertices[61132]={{F(3210778655u),F(1052496443u),F(3199004992u)},{F(3208900786u),F(3193208788u),F(1058796504u)},{F(1044726868u),F(1062489427u),F(1057376315u)},{{F(1033019107u),F(0u)},{F(1024198188u),F(0u)}},{13155,F(1056935908u)},VertexStatus::Valid,0};
 a.vertices[61133]={{F(3210400015u),F(1053014123u),F(3200228514u)},{F(3208959925u),F(3192287929u),F(1058797860u)},{F(1045063154u),F(1062571855u),F(1057212811u)},{{F(1033789577u),F(0u)},{F(1028136468u),F(0u)}},{13143,F(1062677489u)},VertexStatus::Valid,0};
 a.vertices[61261]={{F(3210510398u),F(1053003800u),F(3199724169u)},{F(3209021109u),F(3193516742u),F(1058616521u)},{F(1044535817u),F(1062372063u),F(1057576229u)},{{F(1033937054u),F(0u)},{F(1026276023u),F(0u)}},{13143,F(1060408307u)},VertexStatus::Valid,0};
 a.vertices[61262]={{F(3210130998u),F(1053480497u),F(3200879775u)},{F(3209044897u),F(3192713480u),F(1058654875u)},{F(1044868054u),F(1062459354u),F(1057410709u)},{{F(1034826600u),F(0u)},{F(1031071223u),F(0u)}},{13136,F(1058083680u)},VertexStatus::Valid,0};
 {auto& p=a.patches[40901];p.patch_id=60275;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(958374564u);p.vertex_indices[0]=61132;p.vertex_indices[1]=61133;p.vertex_indices[2]=61261;p.vertex_indices[3]=61262;a.used_patches.push_back(40901);}
 a.vertices[61157]={{F(3209750489u),F(1053880740u),F(1054488788u)},{F(3208996295u),F(3192160952u),F(3206245528u)},{F(3192681836u),F(3210027238u),F(1057245482u)},{{F(1035912021u),F(0u)},{F(3182550056u),F(0u)}},{13139,F(1004386112u)},VertexStatus::Valid,0};
 a.vertices[61286]={{F(3209493661u),F(1054289092u),F(1055006620u)},{F(3208993644u),F(3192855514u),F(3206192388u)},{F(3192488029u),F(3209921292u),F(1057432103u)},{{F(1037266665u),F(0u)},{F(3185331375u),F(0u)}},{13143,F(1052012030u)},VertexStatus::Valid,0};
 a.vertices[61287]={{F(3209858204u),F(1053939408u),F(1054002970u)},{F(3209100479u),F(3193205158u),F(3206024414u)},{F(3192161297u),F(3209831791u),F(1057600218u)},{{F(1035999235u),F(0u)},{F(3180770053u),F(0u)}},{13132,F(1064381046u)},VertexStatus::Valid,0};
 {auto& p=a.patches[40908];p.patch_id=60299;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=61157;p.vertex_indices[1]=61158;p.vertex_indices[2]=61286;p.vertex_indices[3]=61287;a.used_patches.push_back(40908);}
 a.vertices[61263]={{F(3209750489u),F(1053880740u),F(3201972436u)},{F(3208996295u),F(3192160952u),F(1058761880u)},{F(1045198188u),F(1062543590u),F(1057245482u)},{{F(1035912021u),F(0u)},{F(1035066408u),F(0u)}},{13139,F(1004386112u)},VertexStatus::Valid,0};
 a.vertices[61391]={{F(3209858204u),F(1053939408u),F(3201486618u)},{F(3209100479u),F(3193205158u),F(1058540766u)},{F(1044677649u),F(1062348143u),F(1057600218u)},{{F(1035999235u),F(0u)},{F(1033286405u),F(0u)}},{13132,F(1064381046u)},VertexStatus::Valid,0};
 a.vertices[61392]={{F(3209493661u),F(1054289092u),F(3202490268u)},{F(3208993644u),F(3192855514u),F(1058708740u)},{F(1045004381u),F(1062437644u),F(1057432103u)},{{F(1037266665u),F(0u)},{F(1037847727u),F(0u)}},{13143,F(1052012030u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41013];p.patch_id=60404;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=61262;p.vertex_indices[1]=61263;p.vertex_indices[2]=61391;p.vertex_indices[3]=61392;a.used_patches.push_back(41013);}
 a.vertices[61291]={{F(3211287839u),F(1051812111u),F(1049469012u)},{F(3208877836u),F(3196475431u),F(3205948515u)},{F(3190827634u),F(3209441793u),F(1058274580u)},{{F(1032793256u),F(0u)},{F(3164173713u),F(0u)}},{13176,F(1058501128u)},VertexStatus::Valid,0};
 a.vertices[61292]={{F(3211595750u),F(1051131820u),F(1048030344u)},{F(3208737617u),F(3197086811u),F(3205988390u)},{F(3190481712u),F(3209334659u),F(1058445655u)},{{F(1032273218u),F(0u)},{F(3160012205u),F(0u)}},{13197,F(1043534918u)},VertexStatus::Valid,0};
 a.vertices[61420]={{F(3211015156u),F(1052441256u),F(1050283394u)},{F(3209040304u),F(3196516111u),F(3205724985u)},{F(3190661440u),F(3209313598u),F(1058461098u)},{{F(1033629859u),F(0u)},{F(3166745308u),F(0u)}},{13159,F(1049095974u)},VertexStatus::Valid,0};
 a.vertices[61421]={{F(3211335620u),F(1051802207u),F(1049164119u)},{F(3208923070u),F(3197091764u),F(3205745713u)},{F(3190322707u),F(3209201993u),F(1058632626u)},{{F(1033036207u),F(0u)},{F(3163216147u),F(0u)}},{13177,F(1056290256u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41023];p.patch_id=60432;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(967399148u);p.vertex_indices[0]=61291;p.vertex_indices[1]=61292;p.vertex_indices[2]=61420;p.vertex_indices[3]=61421;a.used_patches.push_back(41023);}
 a.vertices[61386]={{F(3211595750u),F(1051131820u),F(3195513992u)},{F(3208737617u),F(3197086811u),F(1058504742u)},{F(1042998064u),F(1061851011u),F(1058445655u)},{{F(1032273218u),F(0u)},{F(1012528557u),F(0u)}},{13197,F(1043534918u)},VertexStatus::Valid,0};
 a.vertices[61387]={{F(3211287839u),F(1051812111u),F(3196952660u)},{F(3208877836u),F(3196475431u),F(1058464867u)},{F(1043343986u),F(1061958145u),F(1058274580u)},{{F(1032793256u),F(0u)},{F(1016690065u),F(0u)}},{13176,F(1058501128u)},VertexStatus::Valid,0};
 a.vertices[61515]={{F(3211335620u),F(1051802207u),F(3196647767u)},{F(3208923070u),F(3197091764u),F(1058262065u)},{F(1042839059u),F(1061718345u),F(1058632626u)},{{F(1033036207u),F(0u)},{F(1015732499u),F(0u)}},{13177,F(1056290256u)},VertexStatus::Valid,0};
 a.vertices[61516]={{F(3211015156u),F(1052441256u),F(3197767042u)},{F(3209040304u),F(3196516111u),F(1058241337u)},{F(1043177792u),F(1061829950u),F(1058461098u)},{{F(1033629859u),F(0u)},{F(1019261660u),F(0u)}},{13159,F(1049095974u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41118];p.patch_id=60527;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(967399148u);p.vertex_indices[0]=61386;p.vertex_indices[1]=61387;p.vertex_indices[2]=61515;p.vertex_indices[3]=61516;a.used_patches.push_back(41118);}
 a.vertices[61414]={{F(3208975500u),F(1054851813u),F(1056188579u)},{F(3208535121u),F(3194305524u),F(3206629972u)},{F(3192621561u),F(3209912201u),F(1057433821u)},{{F(1040470115u),F(0u)},{F(3193465581u),F(0u)}},{13204,F(1055290434u)},VertexStatus::Valid,0};
 a.vertices[61415]={{F(3209249746u),F(1054674479u),F(1055458319u)},{F(3208936307u),F(3193711481u),F(3206191494u)},{F(3192298165u),F(3209817524u),F(1057609701u)},{{F(1038819962u),F(0u)},{F(3188328784u),F(0u)}},{13154,F(1064648835u)},VertexStatus::Valid,0};
 a.vertices[61543]={{F(3208838094u),F(1055102549u),F(1056369031u)},{F(3208154496u),F(3196097531u),F(3206896540u)},{F(3192430252u),F(3209819793u),F(1057594277u)},{{F(1041627375u),F(0u)},{F(3197167713u),F(0u)}},{13257,F(1051238264u)},VertexStatus::Valid,0};
 a.vertices[61544]={{F(3209031196u),F(1055024208u),F(1055826054u)},{F(3208800600u),F(3194808326u),F(3206261476u)},{F(3192111301u),F(3209716858u),F(1057777596u)},{{F(1040400825u),F(0u)},{F(3190798743u),F(0u)}},{13176,F(1064013977u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41126];p.patch_id=60554;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3124931460u);p.vertex_indices[0]=61414;p.vertex_indices[1]=61415;p.vertex_indices[2]=61543;p.vertex_indices[3]=61544;a.used_patches.push_back(41126);}
 a.vertices[61418]={{F(3210316123u),F(1053547844u),F(1052494489u)},{F(3209183059u),F(3195005994u),F(3205744645u)},{F(3191324486u),F(3209524927u),F(1058120531u)},{{F(1035147964u),F(0u)},{F(3174672003u),F(0u)}},{13133,F(1049390816u)},VertexStatus::Valid,0};
 a.vertices[61419]={{F(3210673735u),F(1053024899u),F(1051397868u)},{F(3209129781u),F(3195922228u),F(3205720986u)},{F(3190995015u),F(3209421178u),F(1058290445u)},{{F(1034325151u),F(0u)},{F(3171252620u),F(0u)}},{13144,F(1034647360u)},VertexStatus::Valid,0};
 a.vertices[61547]={{F(3210024771u),F(1054078501u),F(1053150354u)},{F(3209262098u),F(3195340667u),F(3205601851u)},{F(3191150534u),F(3209407869u),F(1058296944u)},{{F(1036311814u),F(0u)},{F(3177594221u),F(0u)}},{13125,F(1064168852u)},VertexStatus::Valid,0};
 a.vertices[61548]={{F(3210381273u),F(1053603940u),F(1052117201u)},{F(3209243343u),F(3196093861u),F(3205542749u)},{F(3190826240u),F(3209299579u),F(1058467942u)},{{F(1035363027u),F(0u)},{F(3173384750u),F(0u)}},{13132,F(1046969428u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41130];p.patch_id=60558;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(957699450u);p.vertex_indices[0]=61418;p.vertex_indices[1]=61419;p.vertex_indices[2]=61547;p.vertex_indices[3]=61548;a.used_patches.push_back(41130);}
 a.vertices[61517]={{F(3210673735u),F(1053024899u),F(3198881516u)},{F(3209129781u),F(3195922228u),F(1058237338u)},{F(1043511367u),F(1061937530u),F(1058290445u)},{{F(1034325151u),F(0u)},{F(1023768972u),F(0u)}},{13144,F(1034647360u)},VertexStatus::Valid,0};
 a.vertices[61518]={{F(3210316123u),F(1053547844u),F(3199978137u)},{F(3209183059u),F(3195005994u),F(1058260997u)},{F(1043840838u),F(1062041279u),F(1058120531u)},{{F(1035147964u),F(0u)},{F(1027188355u),F(0u)}},{13133,F(1049390816u)},VertexStatus::Valid,0};
 a.vertices[61646]={{F(3210381273u),F(1053603940u),F(3199600849u)},{F(3209243343u),F(3196093861u),F(1058059101u)},{F(1043342592u),F(1061815931u),F(1058467942u)},{{F(1035363027u),F(0u)},{F(1025901102u),F(0u)}},{13132,F(1046969428u)},VertexStatus::Valid,0};
 a.vertices[61647]={{F(3210024771u),F(1054078501u),F(3200634002u)},{F(3209262098u),F(3195340667u),F(1058118203u)},{F(1043666886u),F(1061924221u),F(1058296944u)},{{F(1036311814u),F(0u)},{F(1030110573u),F(0u)}},{13125,F(1064168852u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41229];p.patch_id=60657;p.status=PatchCellStatus(1);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(957699450u);p.vertex_indices[0]=61517;p.vertex_indices[1]=61518;p.vertex_indices[2]=61646;p.vertex_indices[3]=61647;a.used_patches.push_back(41229);}
 a.vertices[61521]={{F(3209249746u),F(1054674479u),F(3202941967u)},{F(3208936307u),F(3193711481u),F(1058707846u)},{F(1044814517u),F(1062333876u),F(1057609701u)},{{F(1038819962u),F(0u)},{F(1040845136u),F(0u)}},{13154,F(1064648835u)},VertexStatus::Valid,0};
 a.vertices[61522]={{F(3208975500u),F(1054851813u),F(3203672227u)},{F(3208535121u),F(3194305524u),F(1059146324u)},{F(1045137913u),F(1062428553u),F(1057433821u)},{{F(1040470115u),F(0u)},{F(1045981933u),F(0u)}},{13204,F(1055290434u)},VertexStatus::Valid,0};
 a.vertices[61650]={{F(3209031196u),F(1055024208u),F(3203309702u)},{F(3208800600u),F(3194808326u),F(1058777828u)},{F(1044627653u),F(1062233210u),F(1057777596u)},{{F(1040400825u),F(0u)},{F(1043315095u),F(0u)}},{13176,F(1064013977u)},VertexStatus::Valid,0};
 a.vertices[61651]={{F(3208838094u),F(1055102549u),F(3203852679u)},{F(3208154496u),F(3196097531u),F(1059412892u)},{F(1044946604u),F(1062336145u),F(1057594277u)},{{F(1041627375u),F(0u)},{F(1049684065u),F(0u)}},{13257,F(1051238264u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41233];p.patch_id=60661;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3124931459u);p.vertex_indices[0]=61521;p.vertex_indices[1]=61522;p.vertex_indices[2]=61650;p.vertex_indices[3]=61651;a.used_patches.push_back(41233);}
 a.vertices[61546]={{F(3209667959u),F(1054482164u),F(1054136957u)},{F(3209221319u),F(3194774756u),F(3205715974u)},{F(3191472244u),F(3209512930u),F(1058125722u)},{{F(1037454047u),F(0u)},{F(3181511848u),F(0u)}},{13127,F(1059126757u)},VertexStatus::Valid,0};
 a.vertices[61675]={{F(3209394426u),F(1054935152u),F(1054654605u)},{F(3209211525u),F(3195412682u),F(3205663897u)},{F(3191296657u),F(3209403493u),F(1058291730u)},{{F(1038956349u),F(0u)},{F(3184174906u),F(0u)}},{13132,F(1050287632u)},VertexStatus::Valid,0};
 a.vertices[61676]={{F(3209733209u),F(1054590085u),F(1053751655u)},{F(3209307926u),F(3195766650u),F(3205492306u)},{F(3190979758u),F(3209293376u),F(1058464970u)},{{F(1037628873u),F(0u)},{F(3180287760u),F(0u)}},{13122,F(1062907111u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41236];p.patch_id=60685;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=61546;p.vertex_indices[1]=61547;p.vertex_indices[2]=61675;p.vertex_indices[3]=61676;a.used_patches.push_back(41236);}
 a.vertices[61648]={{F(3209667959u),F(1054482164u),F(3201620605u)},{F(3209221319u),F(3194774756u),F(1058232326u)},{F(1043988596u),F(1062029282u),F(1058125722u)},{{F(1037454047u),F(0u)},{F(1034028200u),F(0u)}},{13127,F(1059126757u)},VertexStatus::Valid,0};
 a.vertices[61776]={{F(3209733209u),F(1054590085u),F(3201235303u)},{F(3209307926u),F(3195766650u),F(1058008658u)},{F(1043496110u),F(1061809728u),F(1058464970u)},{{F(1037628873u),F(0u)},{F(1032804112u),F(0u)}},{13122,F(1062907111u)},VertexStatus::Valid,0};
 a.vertices[61777]={{F(3209394426u),F(1054935152u),F(3202138253u)},{F(3209211525u),F(3195412682u),F(1058180249u)},{F(1043813009u),F(1061919845u),F(1058291730u)},{{F(1038956349u),F(0u)},{F(1036691258u),F(0u)}},{13132,F(1050287632u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41337];p.patch_id=60786;p.status=PatchCellStatus(3);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(2143289344u);p.vertex_indices[0]=61647;p.vertex_indices[1]=61648;p.vertex_indices[2]=61776;p.vertex_indices[3]=61777;a.used_patches.push_back(41337);}
 a.vertices[61803]={{F(3208887365u),F(1055520153u),F(1055764991u)},{F(3208799444u),F(3196400602u),F(3206065711u)},{F(3191432873u),F(3209412796u),F(1058268334u)},{{F(1041395298u),F(0u)},{F(3191759330u),F(0u)}},{13187,F(1061003726u)},VertexStatus::Valid,0};
 a.vertices[61804]={{F(3209140000u),F(1055351119u),F(1055095320u)},{F(3209142687u),F(3196151922u),F(3205668381u)},{F(3191123175u),F(3209297630u),F(1058448594u)},{{F(1040428036u),F(0u)},{F(3187660874u),F(0u)}},{13144,F(1064963120u)},VertexStatus::Valid,0};
 a.vertices[61932]={{F(3208750234u),F(1055771903u),F(1055934188u)},{F(3208435351u),F(3197293688u),F(3206314799u)},{F(3191252576u),F(3209320090u),F(1058408621u)},{{F(1042574578u),F(0u)},{F(3195801560u),F(0u)}},{13238,F(1061038760u)},VertexStatus::Valid,0};
 a.vertices[61933]={{F(3208919869u),F(1055714126u),F(1055439683u)},{F(3208989566u),F(3196708408u),F(3205748493u)},{F(3190950444u),F(3209196437u),F(1058595524u)},{{F(1041403750u),F(0u)},{F(3189999177u),F(0u)}},{13168,F(1064751960u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41446];p.patch_id=60940;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3123443259u);p.vertex_indices[0]=61803;p.vertex_indices[1]=61804;p.vertex_indices[2]=61932;p.vertex_indices[3]=61933;a.used_patches.push_back(41446);}
 a.vertices[61906]={{F(3209140000u),F(1055351119u),F(3202578968u)},{F(3209142687u),F(3196151922u),F(1058184733u)},{F(1043639527u),F(1061813982u),F(1058448594u)},{{F(1040428036u),F(0u)},{F(1040177226u),F(0u)}},{13144,F(1064963120u)},VertexStatus::Valid,0};
 a.vertices[61907]={{F(3208887365u),F(1055520153u),F(3203248639u)},{F(3208799444u),F(3196400602u),F(1058582063u)},{F(1043949225u),F(1061929148u),F(1058268334u)},{{F(1041395298u),F(0u)},{F(1044275682u),F(0u)}},{13187,F(1061003726u)},VertexStatus::Valid,0};
 a.vertices[62035]={{F(3208919869u),F(1055714126u),F(3202923331u)},{F(3208989566u),F(3196708408u),F(1058264845u)},{F(1043466796u),F(1061712789u),F(1058595524u)},{{F(1041403750u),F(0u)},{F(1042515529u),F(0u)}},{13168,F(1064751960u)},VertexStatus::Valid,0};
 a.vertices[62036]={{F(3208750234u),F(1055771903u),F(3203417836u)},{F(3208435351u),F(3197293688u),F(1058831151u)},{F(1043768928u),F(1061836442u),F(1058408621u)},{{F(1042574578u),F(0u)},{F(1048317912u),F(0u)}},{13238,F(1061038760u)},VertexStatus::Valid,0};
 {auto& p=a.patches[41549];p.patch_id=61043;p.status=PatchCellStatus(2);p.incident_area_drop2=F(965708729u);p.signed_solid_angle_sr=F(3123443259u);p.vertex_indices[0]=61906;p.vertex_indices[1]=61907;p.vertex_indices[2]=62035;p.vertex_indices[3]=62036;a.used_patches.push_back(41549);}
 a.offsets.push_back(0);
 a.directions.push_back({F(3201783433u),F(1057325077u),F(3208515119u)});a.original_direction_ids.push_back(11172);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({4369u,6949u,0u,0u,F(1064550614u),F(1063645713u),F(1065352014u),F(855638016u)});
 a.hits.push_back({35124u,54026u,0u,0u,F(1052661621u),F(1054924845u),F(1065347406u),F(0u)});
 a.hits.push_back({35457u,54413u,0u,0u,F(1061633731u),F(1031882991u),F(1065352295u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3201488507u),F(1057730226u),F(3208304832u)});a.original_direction_ids.push_back(11173);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({4370u,6950u,0u,0u,F(1053132394u),F(1018957495u),F(1065350637u),F(855638016u)});
 a.hits.push_back({34910u,53771u,0u,0u,F(1039513476u),F(1064581236u),F(1065351238u),F(864026624u)});
 a.hits.push_back({35237u,54158u,0u,0u,F(1056109542u),F(1062092240u),F(1065351553u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3201488507u),F(1057730226u),F(1060821184u)});a.original_direction_ids.push_back(11236);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({4421u,7001u,0u,0u,F(1058880715u),F(1018957536u),F(1065350636u),F(864026624u)});
 a.hits.push_back({35015u,53876u,0u,0u,F(1063340302u),F(1064581234u),F(1065351238u),F(855638016u)});
 a.hits.push_back({35336u,54257u,0u,0u,F(1057392141u),F(1062092241u),F(1065351551u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3201783433u),F(1057325077u),F(1061031471u)});a.original_direction_ids.push_back(11237);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({4422u,7002u,0u,0u,F(1027863193u),F(1063645714u),F(1065352013u),F(0u)});
 a.hits.push_back({35231u,54133u,0u,0u,F(1059116103u),F(1054924851u),F(1065347404u),F(0u)});
 a.hits.push_back({35558u,54514u,0u,0u,F(1046676727u),F(1031882990u),F(1065352295u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202773326u),F(1057350417u),F(3208203158u)});a.original_direction_ids.push_back(11352);
 a.summaries.push_back({3,3,2,1,0,8,4294967295,0});
 a.hits.push_back({4490u,7078u,0u,0u,F(1062540773u),F(1050298107u),F(1065350107u),F(864026624u)});
 a.hits.push_back({35234u,54155u,0u,0u,F(1035902342u),F(1056097834u),F(1065352684u),F(864026624u)});
 a.hits.push_back({35345u,54284u,0u,0u,F(1007124696u),F(1043917929u),F(1065352853u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202485788u),F(1057745419u),F(3207998138u)});a.original_direction_ids.push_back(11353);
 a.summaries.push_back({3,3,2,1,0,8,4294967295,0});
 a.hits.push_back({4371u,6951u,0u,0u,F(1048647685u),F(1055452079u),F(1065349320u),F(0u)});
 a.hits.push_back({35125u,54027u,0u,0u,F(1062968573u),F(1038023622u),F(1065352596u),F(864026624u)});
 a.hits.push_back({35126u,54028u,0u,0u,F(1059347641u),F(1061628876u),F(1065352608u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3200165883u),F(1060932364u),F(3205534088u)});a.original_direction_ids.push_back(11362);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({3544u,6060u,0u,0u,F(1050791249u),F(1058950679u),F(1065349118u),F(855638016u)});
 a.hits.push_back({33750u,52244u,0u,0u,F(1039869418u),F(1063046597u),F(1065346862u),F(864026624u)});
 a.hits.push_back({34117u,52758u,0u,0u,F(1046377479u),F(1062583234u),F(1065352050u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3200165883u),F(1060932364u),F(1058050440u)});a.original_direction_ids.push_back(11407);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({3583u,6099u,0u,0u,F(1060051287u),F(1058950680u),F(1065349117u),F(855638016u)});
 a.hits.push_back({33837u,52331u,0u,0u,F(1063295810u),F(1063046595u),F(1065346862u),F(864026624u)});
 a.hits.push_back({34200u,52841u,0u,0u,F(1061708541u),F(1062583232u),F(1065352048u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202485788u),F(1057745419u),F(1060514490u)});a.original_direction_ids.push_back(11416);
 a.summaries.push_back({3,3,2,1,0,40,4294967295,0});
 a.hits.push_back({4420u,7000u,0u,0u,F(1061123070u),F(1055452081u),F(1065349320u),F(0u)});
 a.hits.push_back({35229u,54131u,0u,4u,F(1052198542u),F(1061628875u),F(1065352608u),F(864026624u)});
 a.hits.push_back({35230u,54132u,0u,0u,F(1041337354u),F(1038023633u),F(1065352599u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202773326u),F(1057350417u),F(1060719510u)});a.original_direction_ids.push_back(11417);
 a.summaries.push_back({3,3,2,1,0,40,4294967295,0});
 a.hits.push_back({4541u,7129u,0u,0u,F(1043048556u),F(1050298106u),F(1065350106u),F(864026624u)});
 a.hits.push_back({35339u,54260u,0u,4u,F(1063791695u),F(1056097833u),F(1065352684u),F(864026624u)});
 a.hits.push_back({35448u,54387u,0u,0u,F(1065214461u),F(1043917929u),F(1065352853u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3203744820u),F(1057365067u),F(3207876402u)});a.original_direction_ids.push_back(11532);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({4491u,7079u,0u,0u,F(1060636442u),F(1060444159u),F(1065349508u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3203465019u),F(1057749441u),F(3207676899u)});a.original_direction_ids.push_back(11533);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({4372u,6952u,0u,0u,F(1039821555u),F(1063509587u),F(1065351357u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3201207537u),F(1060850631u),F(3205279151u)});a.original_direction_ids.push_back(11542);
 a.summaries.push_back({3,3,2,1,0,40,4294967295,0});
 a.hits.push_back({3662u,6189u,0u,0u,F(1010191417u),F(1050033897u),F(1065351207u),F(864026624u)});
 a.hits.push_back({33930u,52500u,0u,4u,F(1051432012u),F(1052806454u),F(1065352353u),F(855638016u)});
 a.hits.push_back({34022u,52628u,0u,0u,F(1062914247u),F(1052401443u),F(1065352613u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3200391214u),F(1061972045u),F(3203364385u)});a.original_direction_ids.push_back(11546);
 a.summaries.push_back({3,3,3,0,0,32,4294967295,0});
 a.hits.push_back({3315u,5807u,0u,0u,F(1062368697u),F(1063827450u),F(1065351186u),F(0u)});
 a.hits.push_back({33340u,51609u,0u,4u,F(1041459693u),F(1064525701u),F(1065346469u),F(0u)});
 a.hits.push_back({33756u,52250u,0u,0u,F(1060944221u),F(1046460650u),F(1065351797u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3200391214u),F(1061972045u),F(1055880737u)});a.original_direction_ids.push_back(11583);
 a.summaries.push_back({3,3,3,0,0,32,4294967295,0});
 a.hits.push_back({3348u,5840u,0u,0u,F(1043736860u),F(1063827451u),F(1065351186u),F(855638016u)});
 a.hits.push_back({33417u,51686u,0u,4u,F(1062937989u),F(1064525701u),F(1065346469u),F(0u)});
 a.hits.push_back({33831u,52325u,0u,0u,F(1049005380u),F(1046460648u),F(1065351796u),F(855638016u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3201207537u),F(1060850631u),F(1057795503u)});a.original_direction_ids.push_back(11587);
 a.summaries.push_back({3,3,2,1,0,40,4294967295,0});
 a.hits.push_back({3699u,6226u,0u,0u,F(1065166544u),F(1050033898u),F(1065351207u),F(864026624u)});
 a.hits.push_back({34017u,52587u,0u,4u,F(1059730906u),F(1052806454u),F(1065352353u),F(855638016u)});
 a.hits.push_back({34109u,52715u,0u,0u,F(1041554659u),F(1052401447u),F(1065352613u),F(855638016u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3203465019u),F(1057749441u),F(1060193251u)});a.original_direction_ids.push_back(11596);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({4419u,6999u,0u,0u,F(1063301794u),F(1063509586u),F(1065351354u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3203744820u),F(1057365067u),F(1060392754u)});a.original_direction_ids.push_back(11597);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({4540u,7128u,0u,0u,F(1049620940u),F(1060444160u),F(1065349508u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3208664872u),F(3194331746u),F(1058988886u)});a.original_direction_ids.push_back(11628);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({7688u,10453u,0u,0u,F(1049471519u),F(1059736153u),F(1065349002u),F(0u)});
 a.hits.push_back({40787u,60145u,0u,0u,F(1062718712u),F(1051512746u),F(1065352188u),F(0u)});
 a.hits.push_back({41233u,60661u,0u,0u,F(1057916673u),F(1044522755u),F(1065348547u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3208906767u),F(3196524896u),F(1058416164u)});a.original_direction_ids.push_back(11630);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({7804u,10580u,0u,0u,F(1037658194u),F(1065014371u),F(1065352120u),F(864026624u)});
 a.hits.push_back({41118u,60527u,0u,0u,F(1064246554u),F(1047231158u),F(1065352601u),F(864026624u)});
 a.hits.push_back({41549u,61043u,0u,0u,F(1054650286u),F(1052780023u),F(1065348490u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3208906767u),F(3196524896u),F(3205899812u)});a.original_direction_ids.push_back(11679);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({7763u,10539u,0u,0u,F(1063572214u),F(1065014369u),F(1065352120u),F(0u)});
 a.hits.push_back({41023u,60432u,0u,0u,F(1032263484u),F(1047231163u),F(1065352601u),F(855638016u)});
 a.hits.push_back({41446u,60940u,0u,0u,F(1058121769u),F(1052780022u),F(1065348490u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3208664872u),F(3194331746u),F(3206472534u)});a.original_direction_ids.push_back(11681);
 a.summaries.push_back({3,3,3,0,0,0,4294967295,0});
 a.hits.push_back({7645u,10410u,0u,0u,F(1060711151u),F(1059736155u),F(1065349002u),F(864026624u)});
 a.hits.push_back({40688u,60046u,0u,0u,F(1042336799u),F(1051512747u),F(1065352187u),F(864026624u)});
 a.hits.push_back({41126u,60554u,0u,0u,F(1055060484u),F(1044522768u),F(1065348547u),F(847249408u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202232699u),F(1060753944u),F(3205012981u)});a.original_direction_ids.push_back(11722);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({3662u,6189u,0u,0u,F(1060509300u),F(1064960443u),F(1065351195u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3201439944u),F(1061842982u),F(3202911314u)});a.original_direction_ids.push_back(11726);
 a.summaries.push_back({3,3,2,1,0,40,4294967295,0});
 a.hits.push_back({3432u,5936u,0u,0u,F(1054336407u),F(1059981016u),F(1065349072u),F(855638016u)});
 a.hits.push_back({33500u,51865u,0u,4u,F(1035348189u),F(1052452998u),F(1065352830u),F(864026624u)});
 a.hits.push_back({33583u,51993u,0u,0u,F(1058206030u),F(1059846957u),F(1065352556u),F(855638016u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3201439944u),F(1061842982u),F(1055427666u)});a.original_direction_ids.push_back(11763);
 a.summaries.push_back({3,3,2,1,0,40,4294967295,0});
 a.hits.push_back({3463u,5967u,0u,0u,F(1058278708u),F(1059981017u),F(1065349073u),F(855638016u)});
 a.hits.push_back({33577u,51942u,0u,4u,F(1063860964u),F(1052452998u),F(1065352830u),F(864026624u)});
 a.hits.push_back({33660u,52070u,0u,0u,F(1054481764u),F(1059846958u),F(1065352555u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202232699u),F(1060753944u),F(1057529333u)});a.original_direction_ids.push_back(11767);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({3699u,6226u,0u,0u,F(1049875226u),F(1064960445u),F(1065351198u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209008939u),F(3192866531u),F(1058688263u)});a.original_direction_ids.push_back(11808);
 a.summaries.push_back({3,3,2,1,0,8,4294967295,0});
 a.hits.push_back({7687u,10452u,0u,0u,F(1056877062u),F(1041361740u),F(1065349587u),F(864026624u)});
 a.hits.push_back({40901u,60275u,0u,0u,F(1059786324u),F(1060407482u),F(1065352645u),F(872415232u)});
 a.hits.push_back({41013u,60404u,0u,0u,F(1062924264u),F(1064349762u),F(1065352991u),F(847249408u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209243849u),F(3195448181u),F(1058132077u)});a.original_direction_ids.push_back(11810);
 a.summaries.push_back({3,3,2,1,0,8,4294967295,0});
 a.hits.push_back({7803u,10579u,0u,0u,F(1052477713u),F(1053960315u),F(1065348471u),F(847249408u)});
 a.hits.push_back({41229u,60657u,0u,0u,F(1061944971u),F(1062481832u),F(1065352848u),F(864026624u)});
 a.hits.push_back({41337u,60786u,0u,0u,F(1059524058u),F(1062740674u),F(1065352818u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209243849u),F(3195448181u),F(3205615725u)});a.original_direction_ids.push_back(11859);
 a.summaries.push_back({3,3,2,1,0,8,4294967295,0});
 a.hits.push_back({7764u,10540u,0u,0u,F(1059208055u),F(1053960315u),F(1065348469u),F(864026624u)});
 a.hits.push_back({41130u,60558u,0u,0u,F(1045431758u),F(1062481832u),F(1065352848u),F(0u)});
 a.hits.push_back({41236u,60685u,0u,0u,F(1051845706u),F(1062740675u),F(1065352815u),F(0u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209008939u),F(3192866531u),F(3206171911u)});a.original_direction_ids.push_back(11861);
 a.summaries.push_back({3,3,2,1,0,8,4294967295,0});
 a.hits.push_back({7646u,10411u,0u,0u,F(1057008380u),F(1041361742u),F(1065349586u),F(864026624u)});
 a.hits.push_back({40798u,60172u,0u,0u,F(1051321178u),F(1060407479u),F(1065352645u),F(864026624u)});
 a.hits.push_back({40908u,60299u,0u,0u,F(1041514585u),F(1064349764u),F(1065352992u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202471899u),F(1061697756u),F(3202439676u)});a.original_direction_ids.push_back(11906);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({3549u,6065u,0u,0u,F(1019729210u),F(1055413431u),F(1065350628u),F(855638016u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3202471899u),F(1061697756u),F(1054956028u)});a.original_direction_ids.push_back(11943);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({3578u,6094u,0u,0u,F(1064943958u),F(1055413430u),F(1065350628u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209337228u),F(3191384766u),F(1058375321u)});a.original_direction_ids.push_back(11988);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({7568u,10323u,0u,0u,F(1060489608u),F(1058855837u),F(1065348863u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209564870u),F(3193886523u),F(1057836347u)});a.original_direction_ids.push_back(11990);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({7685u,10450u,0u,0u,F(1059045266u),F(1062603218u),F(1065349623u),F(847249408u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209564870u),F(3193886523u),F(3205319995u)});a.original_direction_ids.push_back(12039);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({7648u,10413u,0u,0u,F(1052803293u),F(1062603219u),F(1065349622u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 a.directions.push_back({F(3209337228u),F(3191384766u),F(3205858969u)});a.original_direction_ids.push_back(12041);
 a.summaries.push_back({1,1,1,0,0,0,4294967295,0});
 a.hits.push_back({7529u,10284u,0u,0u,F(1049914609u),F(1058855839u),F(1065348863u),F(864026624u)});
 a.offsets.push_back(a.hits.size());
 return a;
}
}
