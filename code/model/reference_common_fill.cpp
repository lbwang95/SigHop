// All ranks share the same explicit completion model at unobserved edges.
// Model-vs-raw TTQ accuracy without requiring a SigHop index to finish.
// Uses the existing exact CCH, cross-checked against Dijkstra on sampled pairs/worlds.
#include "baselines.hpp"

static void customizeArcs(const Dataset& d, const TreeDecomp& T, const vector<float>& arcs, vector<float>& w) {
    w.assign(2 * T.numShortcuts, FINF);
    for (int64_t e = 0; e < d.m; e++) {
        uint64_t p = T.arcPlace[e];
        if (p != UINT64_MAX) w[p] = min(w[p], arcs[e]);
    }
    for (int v : T.order) {
        int k = (int)T.bag[v].size(); uint64_t base = T.scBase[v], tri = T.triBase[v];
        for (int i = 0; i < k; i++) {
            float toI = w[2*(base+i)+1], fromI = w[2*(base+i)];
            for (int j = i+1; j < k; j++) {
                uint64_t sc = T.triSc[tri++];
                float toJ = w[2*(base+j)+1], fromJ = w[2*(base+j)];
                w[2*sc] = min(w[2*sc], toI+fromJ);
                w[2*sc+1] = min(w[2*sc+1], toJ+fromI);
            }
        }
    }
}
int main(int argc, char** argv) {
    try {
        Args a(argc,argv,1); Dataset d=loadDataset(a.need("data")); RawWorlds raw=loadRaw(a.need("raw"));
        if(raw.N!=d.N)throw runtime_error("model/raw world counts differ");
        Dataset fill=loadDataset(a.need("fill-data"));
        if(fill.n!=d.n || fill.m!=d.m || fill.N!=d.N || fill.from!=d.from || fill.to!=d.to || fill.wday!=d.wday || fill.wslot!=d.wslot || fill.wdow!=d.wdow)throw runtime_error("common fill topology or world tags differ");
        TreeDecomp T;T.build(d); GraphCsr g;g.build(d);CCH cch(d,T);Dijkstra dj(d,g);
        string name=a.get("name",a.need("data"));
        vector<double> alphas;for(const auto& v:splitStr(a.get("alphas","0.5,0.8,0.95"),','))alphas.push_back(stod(v));
        Csv summary(a.need("out"),"dataset,r,method,qset,phi,alpha,reference,mean_rel_err,p95_rel_err,max_rel_err,queries,worlds,mean_raw_coverage,underestimate_rate");
        Csv details(a.need("details"),"dataset,r,qset,s,t,phi,alpha,model_quantile,raw_quantile,rel_err,raw_coverage,worlds");
        Csv validation(a.need("validation"),"dataset,r,qset,phi,checked_model_values,checked_raw_values,mismatches,max_rel_diff");
        for(const auto& ps:splitStr(a.need("phis"),'|')) {
            Phi phi=parsePhi(ps);auto worlds=worldsOf(d,phi);if(worlds.empty())continue;
            int cap=a.geti("max-worlds",0);
            if(cap>0 && worlds.size()>(size_t)cap){mt19937 rng(31);shuffle(worlds.begin(),worlds.end(),rng);worlds.resize(cap);sort(worlds.begin(),worlds.end());}
            vector<size_t> checks(worlds.size());iota(checks.begin(),checks.end(),0);mt19937 rng(997);
            shuffle(checks.begin(),checks.end(),rng);
            int checkWorlds=a.geti("check-worlds",10);
            if(checkWorlds>0)checks.resize(min<size_t>(checkWorlds,checks.size()));
            vector<vector<float>> modelMetrics(worlds.size()), rawMetrics(worlds.size());
            map<size_t,vector<float>> checkRawArcs;
            vector<float> arcs(d.m), builtin;
            for(size_t wi=0;wi<worlds.size();wi++) {
                int w=worlds[wi];
                for(int64_t e=0;e<d.m;e++)arcs[e]=dotf(d.U(e),d.Z(w),(int)d.r);
                customizeArcs(d,T,arcs,modelMetrics[wi]);
                if(wi==0){cch.customize(d.Z(w),builtin);if(builtin!=modelMetrics[wi])throw runtime_error("CCH metric initialization differs");}
                for(int64_t e=0;e<d.m;e++)arcs[e]=dotf(fill.U(e),fill.Z(w),(int)fill.r);
                for(int64_t j=raw.off[w];j<raw.off[w+1];j++)arcs[raw.edge[j]]=raw.w[j];
                customizeArcs(d,T,arcs,rawMetrics[wi]);
                if(find(checks.begin(),checks.end(),wi)!=checks.end())checkRawArcs[wi]=arcs;
            }
            for(const auto& qfile:splitStr(a.need("queries"),',')) {
                string qname=qfile.substr(qfile.find_last_of('/')+1);auto qs=loadQueries(qfile);
                qs.resize(min<size_t>(qs.size(),a.geti("max-queries",200)));if(qs.empty())throw runtime_error("empty accuracy query set");
                vector<size_t> checkedQueries(qs.size());iota(checkedQueries.begin(),checkedQueries.end(),0);
                shuffle(checkedQueries.begin(),checkedQueries.end(),rng);
                int checkQueries=a.geti("check-queries",10);
                if(checkQueries>0)checkedQueries.resize(min<size_t>(checkQueries,checkedQueries.size()));
                int mism=0, checked=0;double maxDiff=0;
                for(size_t qi:checkedQueries)for(size_t wi:checks) {
                    int s=qs[qi].first,t=qs[qi].second,w=worlds[wi];
                    double truth=dj.biQuery(s,t,[&](int e){return dotd(d.U(e),d.Z(w),(int)d.r);});
                    double got=cch.query(s,t,modelMetrics[wi]);
                    double diff=fabs(got-truth)/max(1.0,fabs(truth));maxDiff=max(maxDiff,diff);
                    if(!isfinite(got)||!isfinite(truth)||diff>1e-3)mism++;
                    const auto& weights=checkRawArcs.at(wi);
                    truth=dj.biQuery(s,t,[&](int e){return (double)weights[e];});got=cch.query(s,t,rawMetrics[wi]);
                    diff=fabs(got-truth)/max(1.0,fabs(truth));maxDiff=max(maxDiff,diff);
                    if(!isfinite(got)||!isfinite(truth)||diff>1e-3)mism++;
                    checked++;
                }
                validation.row(joinCsv(name,d.r,qname,phi.name,checked,checked,mism,maxDiff));
                if(mism)throw runtime_error("reference mismatches detected; STOP");
                vector<vector<double>> errors(alphas.size()), coverages(alphas.size());vector<int> underestimated(alphas.size(),0);
                vector<float> mv(worlds.size()),rv(worlds.size());
                for(const auto& q:qs) {
                    for(size_t wi=0;wi<worlds.size();wi++){mv[wi]=cch.query(q.first,q.second,modelMetrics[wi]);rv[wi]=cch.query(q.first,q.second,rawMetrics[wi]);}
                    sort(mv.begin(),mv.end());sort(rv.begin(),rv.end());
                    for(size_t ai=0;ai<alphas.size();ai++) {
                        int rank=quantileRank(alphas[ai],(int)worlds.size());double pred=mv[rank],truth=rv[rank];
                        if(!isfinite(pred)||!isfinite(truth))throw runtime_error("non-finite reference quantile");
                        double error=fabs(pred-truth)/max(1e-9,truth), tolerance=1e-6*max(1.0,pred);
                        double coverage=(double)(upper_bound(rv.begin(),rv.end(),pred+tolerance)-rv.begin())/rv.size();
                        errors[ai].push_back(error);coverages[ai].push_back(coverage);underestimated[ai]+=pred+tolerance<truth;
                        details.row(joinCsv(name,d.r,qname,q.first,q.second,phi.name,alphas[ai],pred,truth,error,coverage,worlds.size()));
                    }
                }
                for(size_t ai=0;ai<alphas.size();ai++) {
                    auto& e=errors[ai];sort(e.begin(),e.end());
                    double mean=accumulate(e.begin(),e.end(),0.0)/e.size();
                    double coverage=accumulate(coverages[ai].begin(),coverages[ai].end(),0.0)/e.size();
                    summary.row(joinCsv(name,d.r,"model-CCH-reference",qname,phi.name,alphas[ai],a.get("reference-name","raw-with-current-model-fill"),
                        mean,e[quantileRank(.95,e.size())],e.back(),e.size(),worlds.size(),coverage,(double)underestimated[ai]/e.size()));
                }
                printf("%s r=%lld %s: queries=%zu worlds=%zu mismatches=%d\n",name.c_str(),(long long)d.r,phi.name.c_str(),qs.size(),worlds.size(),mism);fflush(stdout);
            }
        }
    }catch(const exception& e){fprintf(stderr,"error: %s\n",e.what());return 1;}
    return 0;
}
