// SPDX-License-Identifier: MIT
// Passive independently split quadrature, used by the strict H2 native reader.
using Vec=RANSWallCoordinate;
using Tensor=SolverUtils::RANSTensor;
std::string TensorJSON(const Tensor &v){return J::Rows(v);}
Vec Velocity(const RANSWallEvolutionPatch &p,const Vector &s,const RANSWallGlobalPoint &q,Tensor &g)
{Vec u{};for(int b=0;b<p.Dimension();++b){u[p.PhysicalAxis(b)]=q.Value(s,p.offset[b]);g[p.PhysicalAxis(b)]=q.Gradient(s,p.offset[b]);}return u;}
double Volume(const RANSWallEvolutionPatch &p)
{double result=0;for(const auto &q:p.Quadrature(0,16))result+=q.weight;return result;}
std::string Weak(const RANSWallEvolutionPatch &p,const Vector &s,const Vector &h,double alpha,double time,
    const RANSWallEvolutionData &data)
{
    const auto hash=RANSWallR3::Hash(s);const auto work=p.InitializationWork();
    const double H=p.Parameters().geometry.height;
    double outerSupport=p.upper[2];for(const auto &el:p.fields[0]->Elements())if(el.lower[2]==H)outerSupport=std::min(outerSupport,el.upper[2]);
    std::vector<std::vector<Vector>> tests(3);
    for(int mode=0;mode<3;++mode)for(int block=0;block<p.Gauge();++block)
        tests[mode].push_back(p.fields[block]->Project([&](const Vec &q){
            if(mode==0)return 1.;
            if(mode==1)return q[0]/5;
            return q[2]<H/2 || q[2]>outerSupport?0.:q[2]<H?(q[2]-H/2)/(H/2):(outerSupport-q[2])/(outerSupport-H);}));
    auto contract=[&](const Vector &r)
    {
        std::vector<std::string> modes;
        for(int mode=0;mode<3;++mode)
        {Vector values;for(int b=0;b<p.Gauge();++b){long double v=0;
         for(int j=0;j<p.size[b];++j){v+=(long double)tests[mode][b][j]*r[p.offset[b]+j];}
         values.push_back(v);}
         modes.push_back(J::Numbers(values));}return J::List(modes);
    };
    const auto components=p.ResidualResolutionComponents(s,h,alpha,time,data,false);
    std::vector<std::pair<std::string,std::string>> terms;
    for(const auto &[name,values]:components.terms)terms.push_back({name,contract(values)});
    const auto inner=p.Residual(s,h,alpha,time,data,false,false,0),outer=p.Residual(s,h,alpha,time,data,false,false,1);
    std::vector<std::string> rows;
    auto append=[&](const RANSWallGlobalPoint &v,const RANSWallGlobalPoint &pressure,const RANSWallGlobalPoint &ep,
        double weight,double canonicalWeight,bool in,const EquilibriumAWFProfiles *profile)
    {
        Tensor G{};const auto u=Velocity(p,s,v,G);const auto kgrad=v.Gradient(s,p.offset[p.K()]);
        const auto ge=in?Vec{}:ep.Gradient(s,p.offset[p.Epsilon()]);
        const int pb=in?p.InnerPressure():p.OuterPressure();
        const double k=v.Value(s,p.offset[p.K()]),e=in?profile->Dissipation(v.coordinate[2]):ep.Value(s,p.offset[p.Epsilon()]);
        std::vector<std::string> modes;
        for(int mode=0;mode<3;++mode)
        {const auto tg=v.Gradient(tests[mode][0]);const auto eg=in?Vec{}:ep.Gradient(tests[mode][p.Epsilon()]);
         modes.push_back(J::Object({{"v",J::Number(v.Value(tests[mode][0]))},{"grad_v",J::Numbers(tg)},
          {"p",J::Number(pressure.Value(tests[mode][pb]))},{"e",in?"0":J::Number(ep.Value(tests[mode][p.Epsilon()]))},{"grad_e",J::Numbers(eg)}}));}
        const auto metric=p.Mapping().At(v.coordinate);
        auto atH=v.coordinate;atH[2]=H;Tensor gh{};Vec uh{},gpi{},gkh{};
        if(in){uh=Velocity(p,s,p.Point(0,atH,1),gh);
            gpi=p.Point(p.OuterPressure(),atH,1).Gradient(s,p.offset[p.OuterPressure()]);
            gkh=p.Point(p.K(),atH,1).Gradient(s,p.offset[p.K()]);}
        rows.push_back(J::Object({{"inner",J::Boolean(in)},{"chart",J::Numbers(v.coordinate)},
            {"physical",J::Numbers(p.Physical(v.coordinate))},{"geometry",J::Number(v.geometry)},
            {"element",J::Number(v.element)},{"weight",J::Number(weight)},{"canonical_weight",J::Number(canonicalWeight)},
            {"U",J::Numbers(u)},{"G",TensorJSON(G)},
            {"k",J::Number(k)},{"epsilon",J::Number(e)},{"grad_k",J::Numbers(kgrad)},{"grad_epsilon",J::Numbers(ge)},
            {"pi",J::Number(pressure.Value(s,p.offset[pb]))},{"lambda",J::Number(p.GaugeValue(s))},
            {"history",J::Numbers(Vector{v.Value(h,p.offset[0]),v.Value(h,p.offset[1]),v.Value(h,p.offset[p.K()]),in?0:ep.Value(h,p.offset[p.Epsilon()])})},
            {"K",in?J::Number(profile->K()):"null"},{"cutoff",in?J::Number(profile->Reference().cutoff):"null"},
            {"kappa_eff",in?J::Number(metric.kappa1+metric.kappa2):"null"},
            {"force",J::Numbers(in?p.CartesianForcing(s,v.coordinate,time,data,false):Vec{})},
            {"H_U",in?J::Numbers(uh):"null"},{"H_G",in?TensorJSON(gh):"null"},
            {"H_grad_pi",in?J::Numbers(gpi):"null"},{"H_grad_k",in?J::Numbers(gkh):"null"},{"modes",J::List(modes)}}));
    };
    const auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(24,LibUtilities::eGaussGaussLegendre)];
    for(const auto &el:p.innerVelocity->Elements())if(el.lower[2]==0)
        for(int j=0;j<24;++j)
        {
            Vec anchor{(el.lower[0]+el.upper[0])/2+rule->GetZ()[j]*(el.upper[0]-el.lower[0])/2,0,0};
            const auto selected=p.StateProfile(s,anchor,time,data);auto options=selected->Options();options.quadratureOrder=32;
            const auto metric=p.Mapping().At(anchor);
            auto profile=EquilibriumAWFProfiles::WJ(selected->Geometry(),selected->K(),p.Parameters().model,options,
                {metric.kappa1,metric.kappa2});
            const double wx=rule->GetW()[j]*(el.upper[0]-el.lower[0])/2;
            for(const auto &q:profile.Quadrature())
            {auto a=anchor;a[2]=q.r;append(p.Point(0,a),p.Point(p.InnerPressure(),a),{},
                wx*metric.baseArea*q.weight*q.jacobian,wx*q.weight,true,&profile);}
        }
    for(const auto &v:p.outerVelocity->Quadrature(24))
    {
        const auto a=v.coordinate;
        append(p.Point(0,a),p.Point(p.OuterPressure(),a),p.Point(p.Epsilon(),a),v.weight*p.Mapping().At(a).determinant,v.weight,false,nullptr);
    }
    std::vector<std::string> boundary;
    for(const auto &f:p.Exterior(24))
    {
        std::vector<std::string> modes;
        for(int mode=0;mode<3;++mode){Vector values;
            for(int block=0;block<p.InnerPressure();++block)
                values.push_back(f.inner && block==p.Epsilon()?0:p.GeometryPoint(block,f.face.test.geometry,f.face.test.coordinate).Value(tests[mode][block]));
            modes.push_back(J::Numbers(values));}
        boundary.push_back(J::Object({{"role",J::Quote(f.role)},{"chart",J::Numbers(f.face.test.coordinate)},
            {"geometry",J::Number(f.face.test.geometry)},{"modes",J::List(modes)}}));
    }
    Require(hash==RANSWallR3::Hash(s) && work==p.InitializationWork(),"weak observation mutated state/Lift");
    return J::Object({{"state_sha256",J::Quote(hash)},{"alpha_mass",J::Number(alpha)},
        {"scope",J::Quote("three native polynomial tests:1,x/5,interface_hat; independently split actual-state quadrature; raw loads include prescribed exterior reactions")},
        {"native_inner",contract(inner)},{"native_outer",contract(outer)},
        {"native_terms",J::Object(terms)},{"volume",J::Number(Volume(p))},
        {"geometry_kind",J::Number(p.Mapping().kind==RANSWallNormalMap::Kind::Plane?2:p.Mapping().amplitude==0?0:1)},
        {"samples",J::List(rows)},{"boundary_mode_samples",J::List(boundary)},
        {"governing_solves","0"}});
}
