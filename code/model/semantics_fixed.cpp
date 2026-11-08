#include "baselines.hpp"
int main(int argc,char** argv){
    try{
        Args a(argc,argv,1);Dataset d=loadDataset(a.need("data"));GraphCsr g;g.build(d);Dijkstra dj(d,g);
        if(d.slotSeconds!=900||d.slotsPerDay!=96)throw runtime_error("requires quarter-hour model");
        vector<int> grid(d.days*96,-1);for(int w=0;w<d.N;w++){
            int key=d.wday[w]*96+d.wslot[w];if(grid[key]!=-1)throw runtime_error("duplicate day-slot");grid[key]=w;
        }
        auto at=[&](int day,int slot){return day<0||day>=d.days||slot<0||slot>=96?-1:grid[day*96+slot];};
        string name=a.get("name",a.need("data"));vector<double> alphas{.5,.8,.95};
        Csv out(a.need("out"),"dataset,r,qset,phi,alpha,mean_rel_err,p95_rel_err,max_rel_err,queries,worlds,skipped_queries,reference,status");
        Csv details(a.need("details"),"dataset,r,qset,phi,query_id,s,t,alpha,snapshot_quantile,dynamic_quantile,rel_err,worlds");
        for(string ps:splitStr(a.need("phis"),'|')){
            Phi phi=parsePhi(ps);auto worlds=worldsOf(d,phi);if(worlds.empty())throw runtime_error("NO_DATA");
            for(int w:worlds)if(d.wslot[w]!=d.wslot[worlds[0]])throw runtime_error("mixed departure slots");
            for(string qf:splitStr(a.need("queries"),',')){
                auto qs=loadQueries(qf);qs.resize(min<size_t>(qs.size(),a.geti("max-queries",50)));string qname=qf.substr(qf.find_last_of('/')+1);
                vector<vector<double>> errors(3);int skipped=0;
                for(size_t qi=0;qi<qs.size();qi++){
                    auto [s,t]=qs[qi];vector<double> snap,td;
                    for(int w:worlds){
                        snap.push_back(dj.biQuery(s,t,[&](int e){return dotd(d.U(e),d.Z(w),d.r);}));
                        td.push_back(dj.tdQuery(s,t,d.wday[w],d.wslot[w],at));
                        if(!isfinite(td.back()))break;
                    }
                    if(td.size()!=worlds.size()||!isfinite(td.back())){skipped++;continue;}
                    sort(snap.begin(),snap.end());sort(td.begin(),td.end());
                    for(size_t ai=0;ai<alphas.size();ai++){
                        int rk=quantileRank(alphas[ai],worlds.size());double err=fabs(snap[rk]-td[rk])/td[rk];errors[ai].push_back(err);
                        details.row(joinCsv(name,d.r,qname,phi.name,qi,s,t,alphas[ai],snap[rk],td[rk],err,worlds.size()));
                    }
                }
                for(size_t ai=0;ai<alphas.size();ai++){
                    auto& e=errors[ai];sort(e.begin(),e.end());
                    if(e.empty()){out.row(joinCsv(name,d.r,qname,phi.name,alphas[ai],"NA","NA","NA",0,worlds.size(),skipped,"piecewise_constant_speed","NO_DATA"));continue;}
                    out.row(joinCsv(name,d.r,qname,phi.name,alphas[ai],accumulate(e.begin(),e.end(),0.)/e.size(),e[quantileRank(.95,e.size())],e.back(),e.size(),worlds.size(),skipped,"piecewise_constant_speed","ok"));
                }
            }
        }
    }catch(const exception& e){fprintf(stderr,"error: %s\n",e.what());return 1;}return 0;
}
