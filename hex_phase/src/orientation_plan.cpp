#include <ice_crystal/orientation_plan.hpp>
#include <bit>
#include <charconv>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace iceCrystal
{
std::string sha256(std::string_view data)
{
    constexpr std::uint32_t constants[64]={
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::uint32_t h[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    if(data.size()>(std::numeric_limits<std::uint64_t>::max)()/8-128)throw std::length_error("SHA256 input too long.");
    const auto blocks=(data.size()+9+63)/64;const auto bits=std::uint64_t(data.size())*8;
    for(std::size_t block=0;block<blocks;++block)
    {
        std::uint32_t w[64]{};
        for(unsigned j=0;j<64;++j)
        {
            const auto pos=block*64+j;unsigned byte=0;
            if(pos<data.size())byte=static_cast<unsigned char>(data[pos]);
            else if(pos==data.size())byte=128;
            else if(pos>=blocks*64-8)byte=static_cast<unsigned>((bits>>(8*(blocks*64-1-pos)))&255);
            w[j/4]|=std::uint32_t(byte)<<(24-8*(j%4));
        }
        for(unsigned t=16;t<64;++t)
        {
            const auto a=w[t-15],b=w[t-2];
            w[t]=w[t-16]+(std::rotr(a,7)^std::rotr(a,18)^(a>>3))+w[t-7]+(std::rotr(b,17)^std::rotr(b,19)^(b>>10));
        }
        auto a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],j=h[7];
        for(unsigned t=0;t<64;++t)
        {
            const auto s1=std::rotr(e,6)^std::rotr(e,11)^std::rotr(e,25);
            const auto t1=j+s1+((e&f)^(~e&g))+constants[t]+w[t];
            const auto t2=(std::rotr(a,2)^std::rotr(a,13)^std::rotr(a,22))+((a&b)^(a&c)^(b&c));
            j=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
        }
        h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=j;
    }
    std::string out;out.reserve(64);constexpr char hex[]="0123456789abcdef";
    for(auto v:h)for(int k=7;k>=0;--k)out+=hex[(v>>(4*k))&15];
    return out;
}
namespace
{
std::string_view trim(std::string_view s)
{
    const auto first=s.find_first_not_of(" \t\r");
    if(first==s.npos)return {};
    return s.substr(first,s.find_last_not_of(" \t\r")-first+1);
}
template<class T>T number(std::string_view s)
{
    s=trim(s);T x{};const auto r=std::from_chars(s.data(),s.data()+s.size(),x);
    if(s.empty()||r.ec!=std::errc{}||r.ptr!=s.data()+s.size())throw std::invalid_argument("Malformed orientation number.");
    return x;
}
}
OrientationPlan parse_orientation_plan(std::string_view csv)
{
    OrientationPlan plan;plan.source_sha256=sha256(csv);Sum weights{};bool header=false;
    std::size_t start=0,line_number=0;
    while(start<csv.size())
    {
        const auto end=csv.find('\n',start);const auto length=(end==csv.npos?csv.size():end)-start;
        const auto line=trim(csv.substr(start,length));start=end==csv.npos?csv.size():end+1;++line_number;
        if(line.empty()||line.front()=='#')continue;
        try
        {
            if(!header)
            {
                if(line!="qw,qx,qy,qz,weight,samples,seed")throw std::invalid_argument("Expected CSV header qw,qx,qy,qz,weight,samples,seed.");
                header=true;continue;
            }
            std::array<std::string_view,7> cols;std::size_t pos=0;
            for(unsigned c=0;c<7;++c)
            {
                const auto comma=line.find(',',pos);
                if((c<6&&comma==line.npos)||(c==6&&comma!=line.npos))throw std::invalid_argument("Expected seven CSV fields.");
                cols[c]=line.substr(pos,(comma==line.npos?line.size():comma)-pos);pos=comma+1;
            }
            OrientationNode n;
            for(unsigned c=0;c<4;++c)n.quaternion[c]=number<double>(cols[c]);
            n.rotation=rotation_from_quaternion(n.quaternion[0],n.quaternion[1],n.quaternion[2],n.quaternion[3]);
            n.weight=number<double>(cols[4]);n.samples=number<std::uint64_t>(cols[5]);n.seed=number<std::uint64_t>(cols[6]);
            if(!finite_value(n.weight)||!(n.weight>0)||!n.samples)throw std::invalid_argument("Weights and sample counts must be positive.");
            if(n.samples>(std::numeric_limits<std::uint64_t>::max)()-plan.total_samples)throw std::overflow_error("Planned sample count overflows.");
            plan.total_samples+=n.samples;weights.add(n.weight);plan.nodes.push_back(n);
        }
        catch(const std::exception& e){throw std::invalid_argument("Orientation CSV line "+std::to_string(line_number)+": "+e.what());}
    }
    if(!header||plan.nodes.empty()||!finite_value(weights.value)||::fabs(weights.value-1)>2e-12)
        throw std::invalid_argument("Orientation number-probability weights must sum to one (within 2e-12).");
    plan.input_weight_sum=weights.value;
    // Only remove roundoff in an already-normalized discrete probability measure.
    // Raw angular densities / arbitrary importance weights are NOT silently accepted.
    for(auto& n:plan.nodes)n.weight/=weights.value;
    return plan;
}
OrientationPlan read_orientation_plan(const std::filesystem::path& p)
{
    constexpr std::uintmax_t maximum_bytes=256ull*1024*1024;
    const auto bytes=std::filesystem::file_size(p);
    if(bytes>maximum_bytes)throw std::length_error("Orientation plan exceeds explicit 256 MiB input limit.");
    std::ifstream in(p,std::ios::binary);if(!in)throw std::runtime_error("Cannot read orientation plan.");
    std::string data(static_cast<std::size_t>(bytes),'\0');in.read(data.data(),static_cast<std::streamsize>(bytes));
    if(in.gcount()!=static_cast<std::streamsize>(bytes)||in.bad()||in.peek()!=std::char_traits<char>::eof())
        throw std::runtime_error("Short/changed orientation file.");
    return parse_orientation_plan(data);
}
} // namespace iceCrystal
