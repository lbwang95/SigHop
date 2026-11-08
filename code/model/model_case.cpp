// Exact model-world case profiles with explicit directed paths, independent of SigHop.
#include "baselines.hpp"
#include <set>

int main(int argc, char** argv) {
    try {
        Args a(argc,argv,1); Dataset d=loadDataset(a.need("data"));
        GraphCsr g; g.build(d); Dijkstra dj(d,g); TreeDecomp tree; tree.build(d); CCH cch(d,tree);
        auto qs=loadQueries(a.need("queries")); qs.resize(min<size_t>(qs.size(),a.geti("max-queries",5)));
        if(qs.empty())throw runtime_error("no case queries");
        auto worlds=worldsOf(d,parsePhi(a.get("psi","weekday:dow=0-4")));
        if(worlds.empty())throw runtime_error("no case worlds");
        string name=a.get("name",a.need("data"));
        Csv profile(a.need("out"),"dataset,s,t,slot,worlds,p50,p80,p95,lottr,distinct_selected_paths,slot_seconds");
        Csv perWorld(a.need("details"),"dataset,r,query_id,s,t,world,day,slot,path_id,distance,path_sum,cch_distance,rel_diff");
        Csv paths(a.need("paths"),"dataset,query_id,s,t,path_id,edge_count,edge_ids");
        Csv validation(a.need("validation"),"dataset,r,checked_world_queries,mismatches,max_rel_diff");
        vector<map<vector<int>,int>> ids(qs.size());
        vector<map<int,vector<double>>> times(qs.size());
        vector<map<int,set<int>>> selected(qs.size());
        vector<float> metric; vector<int> arcs;
        long long checked=0,mism=0; double maxDiff=0;
        for(int w:worlds) {
            cch.customize(d.Z(w),metric);
            for(size_t qi=0;qi<qs.size();qi++) {
                auto [s,t]=qs[qi]; auto weight=[&](int e){return dotd(d.U(e),d.Z(w),(int)d.r);};
                double value=dj.pathQuery(s,t,weight,arcs);
                reverse(arcs.begin(),arcs.end()); // pathQuery follows parent arcs backwards.
                int vertex=s; double sum=0; bool valid=isfinite(value);
                for(int e:arcs) {
                    if(e<0 || e>=d.m || d.from[e]!=vertex){valid=false;break;}
                    vertex=d.to[e]; sum+=weight(e);
                }
                valid=valid && vertex==t;
                double ref=cch.query(s,t,metric);
                double diff=max(fabs(sum-value),fabs(ref-value))/max(1.0,fabs(value));
                maxDiff=max(maxDiff,diff); checked++;
                if(!valid || !isfinite(ref) || diff>1e-3) {
                    mism++; validation.row(joinCsv(name,d.r,checked,mism,maxDiff));
                    throw runtime_error("case path continuity or distance mismatch; STOP");
                }
                auto inserted=ids[qi].emplace(arcs,(int)ids[qi].size()); int id=inserted.first->second;
                if(inserted.second) {
                    stringstream edges;
                    for(size_t i=0;i<arcs.size();i++)edges<<(i?";":"")<<arcs[i];
                    paths.row(joinCsv(name,qi,s,t,id,arcs.size(),edges.str()));
                }
                perWorld.row(joinCsv(name,d.r,qi,s,t,w,d.wday[w],d.wslot[w],id,value,sum,ref,diff));
                times[qi][d.wslot[w]].push_back(value); selected[qi][d.wslot[w]].insert(id);
            }
        }
        for(size_t qi=0;qi<qs.size();qi++)for(auto& kv:times[qi]) {
            auto& values=kv.second; sort(values.begin(),values.end()); int k=values.size();
            double p50=values[quantileRank(.5,k)],p80=values[quantileRank(.8,k)],p95=values[quantileRank(.95,k)];
            profile.row(joinCsv(name,qs[qi].first,qs[qi].second,kv.first,k,p50,p80,p95,
                p50>0?p80/p50:numeric_limits<double>::quiet_NaN(),selected[qi][kv.first].size(),d.slotSeconds));
        }
        validation.row(joinCsv(name,d.r,checked,mism,maxDiff));
    } catch(const exception& e) {fprintf(stderr,"error: %s\n",e.what());return 1;}
    return 0;
}
