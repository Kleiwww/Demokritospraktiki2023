//this will accumulate the reco rates from the event and classify them per run
std::map<int,float> reco_events_to_run; 
std::map<int,float> c_reco_events_to_run;
std::map<int,float> first_reco_events_to_run;
std::map<int,float> up_reco_events_to_run;
std::map<int,float> c_up_reco_events_to_run;
std::map<int,float> first_up_reco_events_to_run;
std::map<int,float> evts_3d_muon_to_run;std::map<int,float> evts_3d_shower_to_run;std::map<int,float> evts_mx_shower_to_run; 
std::map<int,float> c_evts_3d_muon_to_run;std::map<int,float> c_evts_3d_shower_to_run;std::map<int,float> c_evts_mx_shower_to_run;
std::map<int,float> first_evts_3d_muon_to_run;std::map<int,float> first_evts_3d_shower_to_run;std::map<int,float> first_evts_mx_shower_to_run;
std::map<int,float> up_evts_3d_muon_to_run;std::map<int,float> up_evts_3d_shower_to_run;std::map<int,float> up_evts_mx_shower_to_run;
std::map<int,float> c_up_evts_3d_muon_to_run;std::map<int,float> c_up_evts_3d_shower_to_run;std::map<int,float> c_up_evts_mx_shower_to_run;
std::map<int,float> first_up_evts_3d_muon_to_run;std::map<int,float> first_up_evts_3d_shower_to_run;std::map<int,float> first_up_evts_mx_shower_to_run;

float icecube_flux_diff2019(float E);
float icecube_flux_diff2021(float E);

void Analysis_arca6_arca8(string file_name = "" ){

  //Set some names (first remove the preproc tag)
  //Input file is something like blahblah_preproc.root
  std::string tempo = file_name;
  std::string to_remove = "_preproc.root";
  size_t m = tempo.find(to_remove);
  size_t n = tempo.find_first_of("_", m + to_remove.length());
  tempo.erase(m, n - m + 1);
  TString tempoTwo = tempo;

  //!!you may change this to match the extra cut!!
  TString RootOut = "./histo_processed/histos_"+ tempoTwo +".root";
 
  std::cout << "===> Analysis Starts ! <===" << std::endl;
  
  //=== Define variables
  int run_number, run_subNumber, num_triggered_hits,pseudo_runid,pseudo_subRunid,num_cherenkov_hits,num_cascade_hits,num_triggered_lines,num_triggered_doms,num_triggered_pmts,num_cherenkov_lines,num_cherenkov_doms,num_cherenkov_pmts,num_cascade_lines,num_cascade_doms,num_cascade_pmts;
  bool flag_muon_3D,flag_shower_3D,flag_shower_MX;
  float livetime,pseudo_livetime,evt_id,evt_num,evt_overlays,reco_type,neutrino_type,jbeta0,jbeta0_deg,D,D3d,cos_zen,fi,R,mean_tres_it,zenith,logEreco,Ereco,E_nu,logEreco2,logbeta0,jlik,GNhit,Snpe,SnpeT,Slen,Elen,best_trk_pos_x,best_trk_pos_y,best_trk_pos_z,best_trk_dir_x,best_trk_dir_y,best_trk_dir_z,Q1,Q1new,diffangle,delta_zenith,mu_pos_x,mu_pos_y,mu_pos_z,nu_pos_x,nu_pos_y,nu_pos_z,cos_zen_mu,zenith_mu,logE_mu,logE_nu,logEbundle,logEdepos,bjorken_y,cos_zen_nu,zenith_nu,w1,w3,w3atm,w3iceHESE,w3loi,w3random,w3ice,w3antares,Ntrack_all,Nallbehind,Nallbehind2,Nallbehind3,NtrackIT10,NtrackIT30,NtrackIT,NtrackEarly,NtrackLate,TrLengthIT,TrLengthIT_2,TrLengthIT_3,ratio130,ratio430,ratio330,ratio110,ratio410,ratio310,ratio1,ratio2,ratio3,ratio4,ratio5,ratio6,ratiol,ratiol_trig,ratio4_trig,myratio50_muon,myratio30_muon,myratio50_cascmuon,myratio30_cascmuon,myratio50_casc,myratio30_casc,myratio50_casc_over_mu,myratio30_casc_over_mu,myratio50_cascmuon_over_mu, myratio30_cascmuon_over_mu,ratio_closehits_casc_over_mu,ratio_closehits_cascmuon_over_mu,redToT_casc_over_mu,redToT_cascmuon_over_mu,diff_theta,ratio_closehits_muon,ratio_closehits_cascmuon,ratio_closehits_casc,redToT_muon,redToT_cascmuon,redToT_casc,diff_dist_mu,diff_dist_casc_mu,diff_dist_casc,max_lik_down,max_lik_up,costheta_min,costheta_max,num_of_good_sol,downsol,upsol,nhits_casc,nhits_casc_100,sum_ToT_casc,min_dist_casc,max_dist_casc,nhits_mu,nhits_mu_100,sum_ToT_mu,min_dist_mu,max_dist_mu,nhits_casc_mu,nhits_casc_mu_100,sum_ToT_casc_mu,min_dist_casc_mu,max_dist_casc_mu,min_diff_sollik,max_diff_sollik,min_diff_sol,max_diff_sol,min_zen_sol,max_zen_sol,diffangle_shower_nu,diffangle_shower_mu,delta_zenith_shower_mu,delta_zenith_shower_nu,zenith_shower,beta0_shower,beta0_shower_deg,lik_shower,Nhit_shower,best_trk_pos_shower_x,best_trk_pos_shower_y,best_trk_pos_shower_z,best_trk_dir_shower_x,best_trk_dir_shower_y,best_trk_dir_shower_z,Ereco_shower,Ereco_shower_corrected,logEreco_shower,dlik,normdlik,itoverlen,delta_zenith_track_shower,diffangle_track_shower,ToT_border_mu,ToT_border_casc,ToT_border_cascmu,ToT_trig,max_ToT_trig,ToT_IT,ToT_allIT,Nborder_hits,Nborder_cherenkov_hits,Nborder_dtres_hits,Nborder_DOMs,Nhits_upper,Nhits_lower,Nhits_border_upper,Nhits_border_lower,Nhits_cherenkov_upper,Nhits_cherenkov_lower,Nhits_border_cherenkov_upper,Nhits_border_cherenkov_lower,ratio_cherenkov_lines,ratio_cherenkov_doms,ratio_cherenkov_pmts,ratio_cascade_doms,ratio_cascade_pmts,NtrackIT50,NtrackIT50_2,NtrackIT30_2,NtrackIT10_2,NtrackIT50_3,NtrackIT30_3,NtrackIT10_3,ratio_cherenkov_hits,ratio_cascade_hits,ratio_cascade_lines;
  float ratio_upper_nhits, ratio_lower_nhits, ratio_border_upper_nhits, ratio_border_lower_nhits;
  float float_num_triggered_hits,float_num_triggered_lines,float_num_triggered_doms,float_num_triggered_pmts,float_num_cherenkov_hits,float_num_cherenkov_lines,float_num_cherenkov_doms,float_num_cherenkov_pmts,float_num_cascade_hits,float_num_cascade_lines,float_num_cascade_doms,float_num_cascade_pmts;
   
  //=====Define ttree
  TFile *input;
  if (!gSystem->AccessPathName( file_name.c_str() )){
    input = TFile::Open( file_name.c_str() );
  }
  else{
    std::cout << "ERROR: could not open file" << std::endl;
    exit(1);
  }
  
  TTree *dataTree = (TTree*)input->Get("ProcessedEvents");
  TTree *InfoTree = (TTree*)input->Get("QualityParameters");
  
  
  //======Initialize variables
  run_number=run_subNumber=num_triggered_hits=pseudo_runid=pseudo_subRunid=num_cherenkov_hits=num_cascade_hits=num_triggered_lines=num_triggered_doms=num_triggered_pmts=num_cherenkov_lines=num_cherenkov_doms=num_cherenkov_pmts=num_cascade_lines=num_cascade_doms=num_cascade_pmts=-1;
  flag_muon_3D=flag_shower_3D=flag_shower_MX=0;
  livetime=pseudo_livetime=evt_id=evt_num=evt_overlays=reco_type=neutrino_type=jbeta0=jbeta0_deg=D=D3d=cos_zen=fi=R=mean_tres_it=zenith=logEreco=Ereco=E_nu=logEreco2=logbeta0=jlik=GNhit=Snpe=SnpeT=Slen=Elen=best_trk_pos_x=best_trk_pos_y=best_trk_pos_z=best_trk_dir_x=best_trk_dir_y=best_trk_dir_z=Q1=Q1new=diffangle=delta_zenith=mu_pos_x=mu_pos_y=mu_pos_z=nu_pos_x=nu_pos_y=nu_pos_z=cos_zen_mu=zenith_mu=logE_mu=logE_nu=logEbundle=logEdepos=bjorken_y=cos_zen_nu=zenith_nu=w1=w3=w3atm=w3iceHESE=w3loi=w3random=w3ice=w3antares=Ntrack_all=Nallbehind=Nallbehind2=Nallbehind3=NtrackIT10=NtrackIT30=NtrackIT=NtrackEarly=NtrackLate=TrLengthIT=TrLengthIT_2=TrLengthIT_3=ratio130=ratio430=ratio330=ratio110=ratio410=ratio310=ratio1=ratio2=ratio3=ratio4=ratio5=ratio6=ratiol=ratiol_trig=ratio4_trig=myratio50_muon=myratio30_muon=myratio50_cascmuon=myratio30_cascmuon=myratio50_casc=myratio30_casc=diff_theta=ratio_closehits_muon=ratio_closehits_cascmuon=ratio_closehits_casc=redToT_muon=redToT_cascmuon=redToT_casc=myratio50_cascmuon_over_mu=myratio30_cascmuon_over_mu=ratio_closehits_cascmuon_over_mu=redToT_cascmuon_over_mu=myratio50_casc_over_mu=myratio30_casc_over_mu=ratio_closehits_casc_over_mu=redToT_casc_over_mu =diff_dist_mu=diff_dist_casc_mu=diff_dist_casc=max_lik_down=max_lik_up=costheta_min=costheta_max=num_of_good_sol=downsol=upsol=nhits_casc=nhits_casc_100=sum_ToT_casc=min_dist_casc=max_dist_casc=nhits_mu=nhits_mu_100=sum_ToT_mu=min_dist_mu=max_dist_mu=nhits_casc_mu=nhits_casc_mu_100=sum_ToT_casc_mu=min_dist_casc_mu=max_dist_casc_mu=min_diff_sollik=max_diff_sollik=min_diff_sol=max_diff_sol=min_zen_sol=max_zen_sol=diffangle_shower_nu=diffangle_shower_mu=delta_zenith_shower_mu=delta_zenith_shower_nu=zenith_shower=beta0_shower=beta0_shower_deg=lik_shower=Nhit_shower=best_trk_pos_shower_x=best_trk_pos_shower_y=best_trk_pos_shower_z=best_trk_dir_shower_x=best_trk_dir_shower_y=best_trk_dir_shower_z=Ereco_shower=Ereco_shower_corrected=logEreco_shower=dlik=normdlik=itoverlen=delta_zenith_track_shower=diffangle_track_shower=ToT_border_mu=ToT_border_casc=ToT_border_cascmu=ToT_trig=max_ToT_trig=ToT_IT=ToT_allIT=Nborder_hits=Nborder_cherenkov_hits=Nborder_dtres_hits=Nborder_DOMs=Nhits_upper=Nhits_lower=Nhits_border_upper=Nhits_border_lower=Nhits_cherenkov_upper=Nhits_cherenkov_lower=Nhits_border_cherenkov_upper=Nhits_border_cherenkov_lower=ratio_cherenkov_lines=ratio_cherenkov_doms=ratio_cherenkov_pmts=ratio_cascade_doms=ratio_cascade_pmts=NtrackIT50=NtrackIT50_2=NtrackIT30_2=NtrackIT10_2=NtrackIT50_3=NtrackIT30_3=NtrackIT10_3=ratio_cherenkov_hits=ratio_cascade_hits=ratio_cascade_lines=-999.9;
  ratio_upper_nhits=ratio_lower_nhits=ratio_border_upper_nhits=ratio_border_lower_nhits=-999.0;
  float_num_triggered_hits=float_num_triggered_lines=float_num_triggered_doms=float_num_triggered_pmts=float_num_cherenkov_hits=float_num_cherenkov_lines=float_num_cherenkov_doms=float_num_cherenkov_pmts=float_num_cascade_hits=float_num_cascade_lines=float_num_cascade_doms=float_num_cascade_pmts=-999.9;
  
  //======================== Read from the ttree and pass to variables ===================
  InfoTree->SetBranchAddress("livetime",&livetime);//run's livetime (duration)
  InfoTree->SetBranchAddress("run_number",&run_number);//run number (int)
  InfoTree->SetBranchAddress("run_subNumber",&run_subNumber);//secondary run number (int)

  dataTree->SetBranchAddress("pseudo_runid", &pseudo_runid); //(int)
  dataTree->SetBranchAddress("pseudo_subRunid", &pseudo_subRunid); //(int)
  dataTree->SetBranchAddress("pseudo_livetime", &pseudo_livetime);
  dataTree->SetBranchAddress("num_triggered_hits", &num_triggered_hits);
  dataTree->SetBranchAddress("evt_id", &evt_id);
  dataTree->SetBranchAddress("evt_num", &evt_num);
  dataTree->SetBranchAddress("evt_overlays", &evt_overlays);
  dataTree->SetBranchAddress("reco_type", &reco_type);
  dataTree->SetBranchAddress("jbeta0", &jbeta0);
  dataTree->SetBranchAddress("jbeta0_deg", &jbeta0_deg);
  dataTree->SetBranchAddress("D", &D);
  dataTree->SetBranchAddress("D3d", &D3d);
  dataTree->SetBranchAddress("cos_zen", &cos_zen);
  dataTree->SetBranchAddress("fi", &fi);
  dataTree->SetBranchAddress("R", &R);
  dataTree->SetBranchAddress("mean_tres_it", &mean_tres_it);
  dataTree->SetBranchAddress("zenith", &zenith);
  dataTree->SetBranchAddress("logEreco", &logEreco);
  dataTree->SetBranchAddress("Ereco", &Ereco);
  dataTree->SetBranchAddress("E_nu", &E_nu);
  dataTree->SetBranchAddress("logEreco2", &logEreco2);
  dataTree->SetBranchAddress("logbeta0", &logbeta0);
  dataTree->SetBranchAddress("jlik", &jlik);
  dataTree->SetBranchAddress("GNhit", &GNhit);
  dataTree->SetBranchAddress("Snpe", &Snpe);
  dataTree->SetBranchAddress("SnpeT", &SnpeT);
  dataTree->SetBranchAddress("Slen", &Slen); //jgandalf
  dataTree->SetBranchAddress("Elen", &Elen); //jenergy
  dataTree->SetBranchAddress("best_trk_pos_x", &best_trk_pos_x);
  dataTree->SetBranchAddress("best_trk_pos_y", &best_trk_pos_y);
  dataTree->SetBranchAddress("best_trk_pos_z", &best_trk_pos_z);
  dataTree->SetBranchAddress("best_trk_dir_x", &best_trk_dir_x);
  dataTree->SetBranchAddress("best_trk_dir_y", &best_trk_dir_y);
  dataTree->SetBranchAddress("best_trk_dir_z", &best_trk_dir_z);
  dataTree->SetBranchAddress("Q1", &Q1); 
  dataTree->SetBranchAddress("Q1new", &Q1new); //lik/ndof
  dataTree->SetBranchAddress("diffangle", &diffangle);
  dataTree->SetBranchAddress("delta_zenith", &delta_zenith);
  dataTree->SetBranchAddress("mu_pos_x", &mu_pos_x);
  dataTree->SetBranchAddress("mu_pos_y", &mu_pos_y);
  dataTree->SetBranchAddress("mu_pos_z", &mu_pos_z);
  dataTree->SetBranchAddress("nu_pos_x", &nu_pos_x);
  dataTree->SetBranchAddress("nu_pos_y", &nu_pos_y);
  dataTree->SetBranchAddress("nu_pos_z", &nu_pos_z);
  dataTree->SetBranchAddress("cos_zen_mu", &cos_zen_mu);
  dataTree->SetBranchAddress("zenith_mu", &zenith_mu);
  dataTree->SetBranchAddress("logE_mu", &logE_mu);
  dataTree->SetBranchAddress("logE_nu", &logE_nu);
  dataTree->SetBranchAddress("logEbundle", &logEbundle);
  dataTree->SetBranchAddress("logEdepos", &logEdepos);
  dataTree->SetBranchAddress("bjorken_y", &bjorken_y);
  dataTree->SetBranchAddress("cos_zen_nu", &cos_zen_nu);
  dataTree->SetBranchAddress("zenith_nu", &zenith_nu);
  dataTree->SetBranchAddress("w3atm", &w3atm);
  dataTree->SetBranchAddress("w3iceHESE", &w3iceHESE);
  dataTree->SetBranchAddress("w3loi", &w3loi);
  dataTree->SetBranchAddress("w3random", &w3random);
  dataTree->SetBranchAddress("w3ice", &w3ice);
  dataTree->SetBranchAddress("w3antares", &w3antares);
  dataTree->SetBranchAddress("w1", &w1);
  dataTree->SetBranchAddress("neutrino_type", &neutrino_type);
  dataTree->SetBranchAddress("Ntrack_all", &Ntrack_all);
  dataTree->SetBranchAddress("Nallbehind", &Nallbehind);
  dataTree->SetBranchAddress("Nallbehind2", &Nallbehind2);
  dataTree->SetBranchAddress("Nallbehind3", &Nallbehind3);
  dataTree->SetBranchAddress("NtrackIT10", &NtrackIT10); //in time with the track
  dataTree->SetBranchAddress("NtrackIT30", &NtrackIT30);
  dataTree->SetBranchAddress("NtrackIT", &NtrackIT);
  dataTree->SetBranchAddress("NtrackEarly", &NtrackEarly);
  dataTree->SetBranchAddress("NtrackLate", &NtrackLate);
  dataTree->SetBranchAddress("TrLengthIT", &TrLengthIT);
  dataTree->SetBranchAddress("TrLengthIT_2", &TrLengthIT_2);
  dataTree->SetBranchAddress("TrLengthIT_3", &TrLengthIT_3);
  dataTree->SetBranchAddress("ratio130", &ratio130);
  dataTree->SetBranchAddress("ratio430", &ratio430);
  dataTree->SetBranchAddress("ratio330", &ratio330);
  dataTree->SetBranchAddress("ratio110", &ratio110);
  dataTree->SetBranchAddress("ratio410", &ratio410);
  dataTree->SetBranchAddress("ratio310", &ratio310);
  dataTree->SetBranchAddress("ratio2", &ratio2);
  dataTree->SetBranchAddress("ratio3", &ratio3);
  dataTree->SetBranchAddress("ratio4", &ratio4);
  dataTree->SetBranchAddress("ratio5", &ratio5);
  dataTree->SetBranchAddress("ratio6", &ratio6);
  dataTree->SetBranchAddress("ratiol", &ratiol);
  dataTree->SetBranchAddress("ratiol_trig", &ratiol_trig);
  dataTree->SetBranchAddress("ratio4_trig", &ratio4_trig);
  dataTree->SetBranchAddress("myratio50_muon", &myratio50_muon); //first hit of track
  dataTree->SetBranchAddress("myratio30_muon", &myratio30_muon);
  dataTree->SetBranchAddress("myratio50_cascmuon", &myratio50_cascmuon);
  dataTree->SetBranchAddress("myratio30_cascmuon", &myratio30_cascmuon);
  dataTree->SetBranchAddress("myratio50_casc", &myratio50_casc);
  dataTree->SetBranchAddress("myratio30_casc", &myratio30_casc);
  dataTree->SetBranchAddress("myratio50_cascmuon_over_mu",&myratio50_cascmuon_over_mu);
  dataTree->SetBranchAddress("myratio30_cascmuon_over_mu",&myratio30_cascmuon_over_mu);
  dataTree->SetBranchAddress("myratio50_casc_over_mu",&myratio50_casc_over_mu);
  dataTree->SetBranchAddress("myratio30_casc_over_mu",&myratio30_casc_over_mu);
  dataTree->SetBranchAddress("redToT_cascmuon_over_mu", &redToT_cascmuon_over_mu);
  dataTree->SetBranchAddress("redToT_casc_over_mu", &redToT_casc_over_mu);
  dataTree->SetBranchAddress("ratio_closehits_cascmuon_over_mu", &ratio_closehits_cascmuon_over_mu);
  dataTree->SetBranchAddress("ratio_closehits_casc_over_mu", &ratio_closehits_casc_over_mu);
  dataTree->SetBranchAddress("diff_theta", &diff_theta);
  dataTree->SetBranchAddress("ratio_closehits_muon", &ratio_closehits_muon);
  dataTree->SetBranchAddress("ratio_closehits_cascmuon", &ratio_closehits_cascmuon);
  dataTree->SetBranchAddress("ratio_closehits_casc", &ratio_closehits_casc);
  dataTree->SetBranchAddress("redToT_muon", &redToT_muon);
  dataTree->SetBranchAddress("redToT_cascmuon", &redToT_cascmuon);
  dataTree->SetBranchAddress("redToT_casc", &redToT_casc);
  dataTree->SetBranchAddress("diff_dist_mu", &diff_dist_mu);
  dataTree->SetBranchAddress("diff_dist_casc_mu", &diff_dist_casc_mu);
  dataTree->SetBranchAddress("diff_dist_casc", &diff_dist_casc);
  dataTree->SetBranchAddress("max_lik_down", &max_lik_down);
  dataTree->SetBranchAddress("max_lik_up", &max_lik_up);
  dataTree->SetBranchAddress("costheta_min", &costheta_min);
  dataTree->SetBranchAddress("costheta_max", &costheta_max);
  dataTree->SetBranchAddress("num_of_good_sol", &num_of_good_sol);
  dataTree->SetBranchAddress("downsol", &downsol);
  dataTree->SetBranchAddress("upsol", &upsol);
  dataTree->SetBranchAddress("nhits_casc", &nhits_casc);
  dataTree->SetBranchAddress("nhits_casc_100", &nhits_casc_100);
  dataTree->SetBranchAddress("sum_ToT_casc", &sum_ToT_casc);
  dataTree->SetBranchAddress("min_dist_casc", &min_dist_casc);
  dataTree->SetBranchAddress("max_dist_casc", &max_dist_casc);
  dataTree->SetBranchAddress("nhits_mu", &nhits_mu);
  dataTree->SetBranchAddress("nhits_mu_100", &nhits_mu_100);
  dataTree->SetBranchAddress("sum_ToT_mu", &sum_ToT_mu);
  dataTree->SetBranchAddress("min_dist_mu", &min_dist_mu);
  dataTree->SetBranchAddress("max_dist_mu", &max_dist_mu);
  dataTree->SetBranchAddress("nhits_casc_mu", &nhits_casc_mu);
  dataTree->SetBranchAddress("nhits_casc_mu_100", &nhits_casc_mu_100);
  dataTree->SetBranchAddress("sum_ToT_casc_mu", &sum_ToT_casc_mu);
  dataTree->SetBranchAddress("min_dist_casc_mu", &min_dist_casc_mu);
  dataTree->SetBranchAddress("max_dist_casc_mu", &max_dist_casc_mu);
  dataTree->SetBranchAddress("min_diff_sollik", &min_diff_sollik);
  dataTree->SetBranchAddress("max_diff_sollik", &max_diff_sollik);
  dataTree->SetBranchAddress("min_diff_sol", &min_diff_sol);
  dataTree->SetBranchAddress("max_diff_sol", &max_diff_sol);
  dataTree->SetBranchAddress("min_zen_sol", &min_zen_sol);
  dataTree->SetBranchAddress("max_zen_sol", &max_zen_sol);
  dataTree->SetBranchAddress("diffangle_shower_nu", &diffangle_shower_nu);
  dataTree->SetBranchAddress("diffangle_shower_mu", &diffangle_shower_mu);
  dataTree->SetBranchAddress("delta_zenith_shower_mu", &delta_zenith_shower_mu);
  dataTree->SetBranchAddress("delta_zenith_shower_nu", &delta_zenith_shower_nu);
  dataTree->SetBranchAddress("zenith_shower", &zenith_shower);
  dataTree->SetBranchAddress("beta0_shower", &beta0_shower);
  dataTree->SetBranchAddress("beta0_shower_deg", &beta0_shower_deg);
  dataTree->SetBranchAddress("lik_shower", &lik_shower);
  dataTree->SetBranchAddress("Nhit_shower", &Nhit_shower);
  dataTree->SetBranchAddress("best_trk_pos_shower_x", &best_trk_pos_shower_x);
  dataTree->SetBranchAddress("best_trk_pos_shower_y", &best_trk_pos_shower_y);
  dataTree->SetBranchAddress("best_trk_pos_shower_z", &best_trk_pos_shower_z);
  dataTree->SetBranchAddress("best_trk_dir_shower_x", &best_trk_dir_shower_x);
  dataTree->SetBranchAddress("best_trk_dir_shower_y", &best_trk_dir_shower_y);
  dataTree->SetBranchAddress("best_trk_dir_shower_z", &best_trk_dir_shower_z);
  dataTree->SetBranchAddress("Ereco_shower", &Ereco_shower);
  dataTree->SetBranchAddress("Ereco_shower_corrected", &Ereco_shower_corrected);
  dataTree->SetBranchAddress("logEreco_shower", &logEreco_shower);
  dataTree->SetBranchAddress("dlik", &dlik);
  dataTree->SetBranchAddress("normdlik", &normdlik);
  dataTree->SetBranchAddress("itoverlen", &itoverlen);
  dataTree->SetBranchAddress("delta_zenith_track_shower", &delta_zenith_track_shower);
  dataTree->SetBranchAddress("diffangle_track_shower", &diffangle_track_shower);
  dataTree->SetBranchAddress("flag_muon_3D", &flag_muon_3D);
  dataTree->SetBranchAddress("flag_shower_3D", &flag_shower_3D);
  dataTree->SetBranchAddress("flag_shower_MX", &flag_shower_MX);
  dataTree->SetBranchAddress("ToT_border_mu", &ToT_border_mu);
  dataTree->SetBranchAddress("ToT_border_casc", &ToT_border_casc);
  dataTree->SetBranchAddress("ToT_border_cascmu", &ToT_border_cascmu);
  dataTree->SetBranchAddress("ToT_trig", &ToT_trig);
  dataTree->SetBranchAddress("max_ToT_trig", &max_ToT_trig);
  dataTree->SetBranchAddress("ToT_IT", &ToT_IT);
  dataTree->SetBranchAddress("ToT_allIT", &ToT_allIT);
  dataTree->SetBranchAddress("Nborder_hits", &Nborder_hits);
  dataTree->SetBranchAddress("Nborder_cherenkov_hits", &Nborder_cherenkov_hits);
  dataTree->SetBranchAddress("Nborder_dtres_hits", &Nborder_dtres_hits);
  dataTree->SetBranchAddress("Nborder_DOMs", &Nborder_DOMs);
  dataTree->SetBranchAddress("Nhits_upper", &Nhits_upper);
  dataTree->SetBranchAddress("Nhits_lower", &Nhits_lower);
  dataTree->SetBranchAddress("Nhits_border_upper", &Nhits_border_upper);
  dataTree->SetBranchAddress("Nhits_border_lower", &Nhits_border_lower);
  dataTree->SetBranchAddress("Nhits_cherenkov_upper", &Nhits_cherenkov_upper);
  dataTree->SetBranchAddress("Nhits_cherenkov_lower", &Nhits_cherenkov_lower);
  dataTree->SetBranchAddress("Nhits_border_cherenkov_upper", &Nhits_border_cherenkov_upper);
  dataTree->SetBranchAddress("Nhits_border_cherenkov_lower", &Nhits_border_cherenkov_lower);
  dataTree->SetBranchAddress("NtrackIT50", &NtrackIT50);
  dataTree->SetBranchAddress("NtrackIT50_2", &NtrackIT50_2);
  dataTree->SetBranchAddress("NtrackIT30_2", &NtrackIT30_2);
  dataTree->SetBranchAddress("NtrackIT10_2", &NtrackIT10_2);
  dataTree->SetBranchAddress("NtrackIT50_3", &NtrackIT50_3);
  dataTree->SetBranchAddress("NtrackIT30_3", &NtrackIT30_3);
  dataTree->SetBranchAddress("NtrackIT10_3", &NtrackIT10_3); 
  dataTree->SetBranchAddress("num_triggered_hits", &num_triggered_hits); //(int)
  dataTree->SetBranchAddress("num_triggered_lines", &num_triggered_lines);//(int)
  dataTree->SetBranchAddress("num_triggered_doms", &num_triggered_doms); //(int)
  dataTree->SetBranchAddress("num_triggered_pmts", &num_triggered_pmts); //(int)
  dataTree->SetBranchAddress ("num_cherenkov_hits", &num_cherenkov_hits); //(int)
  dataTree->SetBranchAddress("num_cherenkov_lines", &num_cherenkov_lines); //(int)
  dataTree->SetBranchAddress("num_cherenkov_doms", &num_cherenkov_doms); //(int)
  dataTree->SetBranchAddress("num_cherenkov_pmts", &num_cherenkov_pmts); //(int)
  dataTree->SetBranchAddress("num_cascade_hits", &num_cascade_hits); //(int)
  dataTree->SetBranchAddress("num_cascade_lines", &num_cascade_lines);  //(int)
  dataTree->SetBranchAddress("num_cascade_doms", &num_cascade_doms); //(int)
  dataTree->SetBranchAddress("num_cascade_pmts", &num_cascade_pmts); //(int)
  dataTree->SetBranchAddress("ratio_cherenkov_hits", &ratio_cherenkov_hits);
  dataTree->SetBranchAddress("ratio_cherenkov_lines", &ratio_cherenkov_lines);
  dataTree->SetBranchAddress("ratio_cherenkov_doms", &ratio_cherenkov_doms);
  dataTree->SetBranchAddress("ratio_cherenkov_pmts", &ratio_cherenkov_pmts);
  dataTree->SetBranchAddress("ratio_cascade_hits", &ratio_cascade_hits);
  dataTree->SetBranchAddress("ratio_cascade_lines", &ratio_cascade_lines);
  dataTree->SetBranchAddress("ratio_cascade_doms", &ratio_cascade_doms);
  dataTree->SetBranchAddress("ratio_cascade_pmts", &ratio_cascade_pmts);
  
  //========= new variables (defined from the loaded ones) ================//

  //=====================================================================//
 
  //===================== Define plots =============================================================
  TH1D CosZen ("CosZen","; Cos(reco_zenith)",100,-1,1);     
  TH1D Lik ("Lik","; Likelihood", 700,-100,600);
  TH1D LogBeta0 ("LogBeta0","; LogBeta0",400,-4,0);
  TH1D Beta0_deg ("Beta0_deg","; LogBeta0deg",600,-4,2);
  TH1D DOMs ("DOMs","DOMs; # of DOMs",200,0,200);
  TH1D Cherenkov_DOMs ("Cherenkov_DOMs","Cherenkov_DOMs; # of Cherenkov DOMs",200,0,200);
  TH1D Q1value ("Q1value","Q1value; -Likelihood/Nhits",500,-5,5);
  TH1D NNhits ("NNhits","NNhits; # of Nhits",1000,0,1000);
  TH1D LogEreco ("LogEreco","; LogEreco ",100,0,8);
  TH1D TrLen  ("TrLen","; Track Length",1000,0,1000);
  TH1D Npe  ("Npe","; # of p.e",10000,0,10000);
  TH1D Ratiopmts  ("Ratiopmts","; Chere.PMTs/Reco.PMTs",100,0,1); 
  TH1D Ratiodoms  ("Ratiodoms","; Chere.DOMs/Reco.DOMs",100,0,1);
  TH2D Allreco_rate ("Allreco_rate","; run_nr",4667,9333,14000,10000,0,100);
  TH1D Rvalue ("Rvalue","R; R(m)",1000,0,1000);
  TH1D All_triggered_hits ("All_triggered_hits","All_triggered_hits; # of hits",1000,0,1000);
  TH1D Phi ("Phi","Phi; Azimuth [deg]",360,-180,180);
  TH1D Zenith_shower ("Zenith_shower","zenith_shower; [deg]",360,0,360);
  TH1D DLik ("DLik","; DLikelihood", 700,-100,600);
  TH1D Ratio6 ("Ratio6","; Ratio6",1200,0,1.2);
  TH1D MaxZenSol ("MaxZenSol","MaxZenSol; [deg]",360,0,360);
  TH1D RedToT_muon ("RedToT_muon","; Total_tot/nhits",300,0,300);
  TH1D Zenith ("Zenith","Zenith; [deg]",360,0,360);
  TH1D Delta_ZeniTh_track_shower ("Delta_ZeniTh_track_shower","Delta_ZeniTh_track_shower; [deg]",360,-180,180);
  TH1D MinZenSol ("MinZenSol","MinZenSol; [deg]",360,0,360);
  TH1D ItoLen ("ItoLen","; ItoLen",1200,0,1.2);
  TH1D Ratio330 ("Ratio330","; Ratio330",1200,0,1.2);
  TH1D Num_good_sol ("Num_good_sol",";# of good solutions",100,0,100);
  TH1D PosZ  ("PosZ","; PosZ",1000,0,1000);
  TH1D LogEnu ("LogEnu","; LogEnu ",100,0,8);
  TH1D LogAEnergy ("LogAEnergy","; LogAEnergy ",100,0,8);
  TH1D NLines ("NLines","; # of lines",12,0,12);
  TH1D LogEreco_cor ("LogEreco_cor","; LogEreco_cor ",100,0,8);
  TH1D Chere_NLines ("Chere_NLines","; # of lines",12,0,12);
  //TH1D Reco_Dclosest ("Reco_Dclosest","; d_closest ", 1000, 0, 1000);
  //TH1D Reco_tres ("Reco_tres","; reco_tres ", 1000, -500, 500);
  //TH1D Reco_Cos_Angle ("Reco_Cos_Angle","; reco_cos_angle ", 500, -360, 360);
  TH1D Ratio_upper_hemisphere ("Ratio_upper_hemisphere","; upper_hits/all_hits",100,0,1);
  TH1D Ratio_lower_hemisphere ("Ratio_lower_hemisphere","; lowerr_hits/all_hits",100,0,1);
  TH1D Ratio_border_upper_hemisphere ("Ratio_border_upper_hemisphere","; border_upper_hits/all_border_hits",100,0,1);
  TH1D Ratio_border_lower_hemisphere ("Ratio_border_lower_hemisphere","; border_lowerr_hits/all_border_hits",100,0,1);
  TH1D Ratio1 ("Ratio1","; Ratio1",1200,0,1.2);
  TH1D Ratio110 ("Ratio110","; Ratio110",1200,0,1.2);
  TH1D Ratio130 ("Ratio130","; Ratio130",1200,0,1.2);
  TH1D NtrackIT_ratio ("NtrackIT_ratio","; NtrackIT/Nhits",1200,0,1.2);
  TH1D TrLenIT_3  ("TrLenIT_3","; Track Length",1000,0,1000);

  TH1D  h_ToT_border_mu("h_ToT_border_mu",";ToT_border_mu ",400,0,4000);
  TH1D  h_ToT_border_casc("h_ToT_border_casc",";ToT_border_casc ",400,0,4000);
  TH1D  h_ToT_border_cascmu("h_ToT_border_cascmu",";ToT_border_cascmu ",400,0,4000);
  TH1D  h_ToT_trig("h_ToT_trig",";ToT_trig ",400,0,4000);
  TH1D  h_max_ToT_trig("h_max_ToT_trig",";max_ToT_trig ",400,0,4000);
  TH1D  h_ToT_IT("h_ToT_IT",";ToT_IT ",400,0,4000);
  TH1D  h_ToT_allIT("h_ToT_allIT",";ToT_allIT ",400,0,4000);
  TH1D  h_sum_ToT_casc("h_sum_ToT_casc",";sum_ToT_casc" ,400,0,4000);
  TH1D  h_sum_ToT_mu("h_sum_ToT_mu",";sum_ToT_mu" ,400,0,4000);
  TH1D  h_sum_ToT_casc_mu("h_sum_ToT_casc_mu",";sum_ToT_casc_mu" ,400,0,4000);
  //----------------------------------------------------------
  TH1D h_Ratio430 ("h_Ratio430","; Ratio430",1200,0,1.2);
  TH1D h_Ratio5 ("h_Ratio5","; Ratio5",1200,0,1.2);
  TH1D h_DiffDistance_mu  ("h_DiffDistance_mu","; dmax-dmin [m]",1400,0,1400);
  TH1D h_RatioCloseHits_mu  ("h_RatioCloseHits_mu","; nhits_400m/nhits",400,0,1);
  TH1D h_NtrackEarly("h_NtrackEarly", ";NtrackEarly",400,0,4000);
  TH1D h_D("h_D", ";D", 400, 0, 4000);
  TH1D h_NtrackIT30("h_NtrackIT30", ";NtrackIT30", 400, 0, 4000);
  TH1D h_NtrackLate("h_NtrackLate", ";NtrackLate", 400, 0, 4000);
  TH1D h_myratio50_muon("h_myratio50_muon", ";myratio50_muon", 400, 0, 1);
  TH1D h_ratio_cherenkov_lines("h_ratio_cherenkov_lines", ";ratio_cherenkov_lines", 400, 0, 1);
  TH1D h_SnpeT("h_SnpeT", ";SnpeT", 400, 0, 4000); //!!Same with Npe -- to be deleted
  TH1D h_mean_tres_it("h_mean_tres_it", ";mean_tres_it", 400, 0, 4000);
  TH1D h_max_lik_up("h_max_lik_up", ";max_lik_up", 400, 0, 4000);
  TH1D h_max_lik_down("h_max_lik_down", ";max_lik_down", 400, 0, 4000);
  TH1D h_diff_theta("h_diff_theta", ";diff_theta", 400, 0, 4000);
  TH1D h_diff_dist_casc_mu("h_diff_dist_casc_mu", ";diff_dist_casc_mu", 400, 0, 4000);
  TH1D h_diff_dist_casc("h_diff_dist_casc", ";diff_dist_casc", 400, 0, 4000);
  TH1D h_ratio_closehits_cascmuon("h_ratio_closehits_cascmuon", ";ratio_closehits_cascmuon", 400, 0, 1);
  TH1D h_ratio_closehits_casc("h_ratio_closehits_casc", ";ratio_closehits_casc", 400, 0,1);
  TH1D h_redToT_cascmuon("h_redToT_cascmuon", ";redToT_cascmuon", 400, 0, 4000);
  TH1D h_redToT_casc("h_redToT_casc", ";redToT_casc", 400, 0, 4000);
  TH1D h_myratio50_cascmuon("h_myratio50_cascmuon", ";myratio50_cascmuon", 400, 0, 1);
  TH1D h_myratio50_casc("h_myratio50_casc", ";myratio50_casc", 400, 0, 1);
  TH1D h_myratio30_cascmuon("h_myratio30_cascmuon", ";myratio30_cascmuon", 400, 0, 1);
  TH1D h_myratio30_muon("h_myratio30_muon", ";myratio30_muon", 400, 0, 1);
  TH1D h_myratio30_casc("h_myratio30_casc", ";myratio30_casc", 400, 0, 1);
  TH1D h_min_diff_sollik("h_min_diff_sollik",";min_diff_sollik", 400, 0, 4000);
  TH1D h_beta0_shower_deg("h_beta0_shower_deg", ";beta0_shower_deg", 400, 0, 4000);
  TH1D h_lik_shower("h_lik_shower", ";lik_shower", 400, 0, 4000);
  TH1D h_best_trk_pos_shower_z("h_best_trk_pos_shower_z", ";best_trk_pos_shower_z", 400, 0, 4000);
  TH1D h_normdlik("h_normdlik", ";normdlik", 400, 0, 4000);
  TH1D h_upsol("h_upsol", ";upsol", 100, 0, 60);
  TH1D h_Slen("h_Slen", ";Slen", 1000, 0, 1000);
  TH1D h_max_diff_sollik("h_max_diff_sollik", ";max_diff_sollik", 400, 0, 4000);
  TH1D h_diffangle_track_shower("h_diffangle_track_shower", ";diffangle_track_shower", 400, 0, 4000);
  TH1D h_myratio50_cascmuon_over_mu("h_myratio50_cascmuon_over_mu", ";myratio50_cascmuon_over_mu",400, 0, 1);
  TH1D h_myratio30_cascmuon_over_mu("h_myratio30_cascmuon_over_mu", ";myratio30_cascmuon_over_mu",400, 0, 1);
  TH1D h_myratio50_casc_over_mu("h_myratio50_casc_over_mu", ";myratio50_casc_over_mu",400, 0, 1);
  TH1D h_myratio30_casc_over_mu("h_myratio30_casc_over_mu", ";myratio30_casc_over_mu",400, 0, 1);
  TH1D h_ratio_closehits_cascmuon_over_mu("h_ratio_closehits_cascmuon_over_mu", ";ratio_closehits_cascmuon_over_mu",400, 0, 1);
  TH1D h_ratio_closehits_casc_over_mu("h_ratio_closehits_casc_over_mu", ";ratio_closehits_casc_over_mu",400, 0, 1);
  TH1D h_redToT_cascmuon_over_mu("h_redToT_cascmuon_over_mu", ";redToT_cascmuon_over_mu",400, 0, 1);
  TH1D h_redToT_casc_over_mu("h_redToT_casc_over_mu", ";redToT_casc_over_mu",400, 0, 1);
  TH1D h_bjorken_y("h_bjorken_y", ";bjorken_y",400, 0, 1);
  TH1D h_diffangle("h_diffangle", ";diffangle",180, 0, 180);
  TH1D h_LogE_mu ("h_LogE_mu","; LogE_mu ",100,0,10);
  TH1D h_LogEbundle ("h_LogEbundle","; LogEbundle ",100,0,10);
  //TH1Dh_logE_mu_max ("//h_logE_mu_max","; LogE_mu_max ",100,0,10);
  TH1D h_cos_zen_mu ("h_cos_zen_mu","; Cos(zenith_mu)",100,-1,1);
  //--------------------------------------------------------- 
 
  TH2D Lik_vs_NNhits ("Lik_vs_NNhits","; ",1000,0,1000,700,-100,600);
  TH2D Lik_vs_Q1value ("Lik_vs_Q1value","; ",500,-5,5,700,-100,600);
  TH2D NNhits_vs_Q1value ("NNhits_vs_Q1value","; ",500,-5,5,1000,0,1000);
  TH2D Nlines_vs_zenith ("Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D Chere_Nlines_vs_zenith ("Chere_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D TrLengthIT_2_vs_Slen ("TrLengthIT_2_vs_Slen",";Slen;TrLengthIT_2",500,0,1000,500,0,1000);
  TH2D Diffangle_vs_run ("Diffangle_vs_run",";diffangle;run",4667,9333,14000,150,0,150);
  TH2D LogEreco_vs_CosZen ("LogEreco_vs_CosZen",";CosZen;LogEreco",100,-1,1,100,0,9);
  TH2D LogEreco_vs_LogEreco2 ("LogEreco_vs_LogEreco2",";LogEreco2;LogEreco",100,-1,10,100,-1,10);
  TH2D LogEreco_vs_LogE_mu ("LogEreco_vs_LogE_mu",";LogE_mu;LogEreco",100,-1,10,100,-1,10);
  TH2D LogEreco2_vs_LogE_mu ("LogEreco2_vs_LogE_mu",";LogE_mu;LogEreco2",100,-1,10,100,-1,10);
  TH2D Rate3Dmuon ("Rate3Dmuon","; run_nr",4667,9333,14000,10000,0,100);
  TH2D Rate3Dshower ("Rate3Dshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D RateMXshower ("RateMXshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D LogEreco_vs_LogE_neu ("LogEreco_vs_LogE_neu",";LogE_neu;LogEreco",100,-1,10,100,-1,10);
  TH2D LogEreco2_vs_LogE_neu ("LogEreco2_vs_LogE_neu",";LogE_neu;LogEreco2",100,-1,10,100,-1,10);
  TH2D LogEresolution_mu_vs_CosZen ("LogEresolution_mu_vs_CosZen",";CosZen;(LogEreco-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D LogEresolution_mu_cor_vs_CosZen ("LogEresolution_mu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D LogEresolution_neu_vs_CosZen ("LogEresolution_neu_vs_CosZen",";CosZen;(LogEreco-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  TH2D LogEresolution_neu_cor_vs_CosZen ("LogEresolution_neu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  
  //Anti-noise cuts
  TH1D C_CosZen ("C_CosZen","C_CosZen; Cos(Zenith)",100,-1,1);
  TH1D C_Lik ("C_Lik","C_Lik; -Lik",700,-100,600);
  TH1D C_LogBeta0 ("C_LogBeta0","C_LogBeta0; Log(Beta0)",400,-4,0);
  TH1D C_Beta0_deg ("C_Beta0_deg","; LogBeta0deg",600,-4,2);
  TH1D C_DOMs ("C_DOMs","DOMs; # of DOMs",200,0,200);
  TH1D C_Cherenkov_DOMs ("C_Cherenkov_DOMs","Cherenkov_DOMs; # of Cherenkov DOMs",200,0,200);
  TH1D C_Q1value ("C_Q1value","C_Q1value; -Likelihood/Nhits",500,-5,5);
  TH1D C_NNhits ("C_NNhits","C_NNhits; # of Nhits",1000,0,1000);
  TH1D C_LogEreco ("C_LogEreco","; LogEreco ",100,0,8);
  TH1D C_TrLen  ("C_TrLen","; Track Length",1000,0,1000);
  TH1D C_Npe  ("C_Npe","; # of p.e",10000,0,10000);
  TH1D C_Ratiopmts  ("C_Ratiopmts","; Chere.PMTs/Reco.PMTs",100,0,1); 
  TH1D C_Ratiodoms  ("C_Ratiodoms","; Chere.DOMs/Reco.DOMs",100,0,1);
  TH2D C_Allreco_rate ("C_Allreco_rate","; run_nr",4667,9333,14000,10000,0,100);
  TH1D C_Rvalue ("C_Rvalue","C_R; R(m)",1000,0,1000);
  TH1D C_All_triggered_hits ("C_All_triggered_hits","C_All_triggered_hits; # of hits",1000,0,1000);
  TH1D C_Phi ("C_Phi","C_Phi; Azimuth [deg]",360,-180,180);
  TH1D C_Zenith_shower ("C_Zenith_shower","C_zenith_shower; [deg]",360,0,360);
  TH1D C_DLik ("C_DLik","; C_DLikelihood", 700,-100,600);
  TH1D C_Ratio6 ("C_Ratio6","; C_Ratio6",1200,0,1.2);
  TH1D C_MaxZenSol ("C_MaxZenSol","C_MaxZenSol; [deg]",360,0,360);
  TH1D C_RedToT_muon ("C_RedToT_muon","; Total_tot/nhits",300,0,300);
  TH1D C_Zenith ("C_Zenith","C_Zenith; [deg]",360,0,360);
  TH1D C_Delta_ZeniTh_track_shower ("C_Delta_ZeniTh_track_shower","C_Delta_ZeniTh_track_shower; [deg]",360,-180,180);
  TH1D C_MinZenSol ("C_MinZenSol","C_MinZenSol; [deg]",360,0,360);
  TH1D C_ItoLen ("C_ItoLen","; C_ItoLen",1200,0,1.2);
  TH1D C_Ratio330 ("C_Ratio330","; C_Ratio330",1200,0,1.2);
  TH1D C_Num_good_sol ("C_Num_good_sol",";# of good solutions",100,0,100);
  TH1D C_PosZ  ("C_PosZ","; C_PosZ",1000,0,1000);
  TH1D C_LogEnu ("C_LogEnu","; C_LogEnu ",100,0,8);
  TH1D C_LogAEnergy ("C_LogAEnergy","; C_LogAEnergy ",100,0,8);
  TH1D C_NLines ("C_NLines","; # of lines",12,0,12);
  TH1D C_LogEreco_cor ("C_LogEreco_cor","; LogEreco_cor ",100,0,8);
  TH1D C_Chere_NLines ("C_Chere_NLines","; # of lines",12,0,12);
  //TH1D C_Reco_Dclosest ("C_Reco_Dclosest","; d_closest ", 1000, 0, 1000);
  //TH1D C_Reco_tres ("C_Reco_tres","; reco_tres ", 1000, -500, 500);
  //TH1D C_Reco_Cos_Angle ("C_Reco_Cos_Angle","; reco_cos_angle ", 500, -360, 360);
  TH1D C_Ratio_upper_hemisphere ("C_Ratio_upper_hemisphere","; upper_hits/all_hits",100,0,1);
  TH1D C_Ratio_lower_hemisphere ("C_Ratio_lower_hemisphere","; lowerr_hits/all_hits",100,0,1);
  TH1D C_Ratio_border_upper_hemisphere ("C_Ratio_border_upper_hemisphere","; border_upper_hits/all_border_hits",100,0,1);
  TH1D C_Ratio_border_lower_hemisphere ("C_Ratio_border_lower_hemisphere","; border_lowerr_hits/all_border_hits",100,0,1);
  TH1D C_Ratio1 ("C_Ratio1","; Ratio1",1200,0,1.2);
  TH1D C_Ratio110 ("C_Ratio110","; Ratio110",1200,0,1.2);
  TH1D C_Ratio130 ("C_Ratio130","; Ratio130",1200,0,1.2);
  TH1D C_NtrackIT_ratio ("C_NtrackIT_ratio","; NtrackIT/Nhits",1200,0,1.2);
  TH1D C_TrLenIT_3  ("C_TrLenIT_3","; Track Length",1000,0,1000);

  TH1D  C_h_ToT_border_mu("C_h_ToT_border_mu",";ToT_border_mu ",400,0,4000);
  TH1D  C_h_ToT_border_casc("C_h_ToT_border_casc",";ToT_border_casc ",400,0,4000);
  TH1D  C_h_ToT_border_cascmu("C_h_ToT_border_cascmu",";ToT_border_cascmu ",400,0,4000);
  TH1D  C_h_ToT_trig("C_h_ToT_trig",";ToT_trig ",400,0,4000);
  TH1D  C_h_max_ToT_trig("C_h_max_ToT_trig",";max_ToT_trig ",400,0,4000);
  TH1D  C_h_ToT_IT("C_h_ToT_IT",";ToT_IT ",400,0,4000);
  TH1D  C_h_ToT_allIT("C_h_ToT_allIT",";ToT_allIT ",400,0,4000);
  TH1D  C_h_sum_ToT_casc("C_h_sum_ToT_casc",";sum_ToT_casc" ,400,0,4000);
  TH1D  C_h_sum_ToT_mu("C_h_sum_ToT_mu",";sum_ToT_mu" ,400,0,4000);
  TH1D  C_h_sum_ToT_casc_mu("C_h_sum_ToT_casc_mu",";sum_ToT_casc_mu" ,400,0,4000);
  //----------------------------------------------------------
  TH1D C_h_Ratio430 ("C_h_Ratio430","; Ratio430",1200,0,1.2);
  TH1D C_h_Ratio5 ("C_h_Ratio5","; Ratio5",1200,0,1.2);
  TH1D C_h_DiffDistance_mu  ("C_h_DiffDistance_mu","; dmax-dmin [m]",1400,0,1400);
  TH1D C_h_RatioCloseHits_mu  ("C_h_RatioCloseHits_mu","; nhits_400m/nhits",400,0,1);
  TH1D C_h_NtrackEarly("C_h_NtrackEarly", ";NtrackEarly",400,0,4000);
  TH1D C_h_D("C_h_D", ";D", 400, 0, 4000);
  TH1D C_h_NtrackIT30("C_h_NtrackIT30", ";NtrackIT30", 400, 0, 4000);
  TH1D C_h_NtrackLate("C_h_NtrackLate", ";NtrackLate", 400, 0, 4000);
  TH1D C_h_myratio50_muon("C_h_myratio50_muon", ";myratio50_muon", 400, 0, 1);
  TH1D C_h_ratio_cherenkov_lines("C_h_ratio_cherenkov_lines", ";ratio_cherenkov_lines", 400, 0, 1);
  TH1D C_h_SnpeT("C_h_SnpeT", ";SnpeT", 400, 0, 4000);
  TH1D C_h_mean_tres_it("C_h_mean_tres_it", ";mean_tres_it", 400, 0, 4000);
  TH1D C_h_max_lik_up("C_h_max_lik_up", ";max_lik_up", 400, 0, 4000);
  TH1D C_h_max_lik_down("C_h_max_lik_down", ";max_lik_down", 400, 0, 4000);
  TH1D C_h_diff_theta("C_h_diff_theta", ";diff_theta", 400, 0, 4000);
  TH1D C_h_diff_dist_casc_mu("C_h_diff_dist_casc_mu", ";diff_dist_casc_mu", 400, 0, 4000);
  TH1D C_h_diff_dist_casc("C_h_diff_dist_casc", ";diff_dist_casc", 400, 0, 4000);
  TH1D C_h_ratio_closehits_cascmuon("C_h_ratio_closehits_cascmuon", ";ratio_closehits_cascmuon", 400, 0, 1);
  TH1D C_h_ratio_closehits_casc("C_h_ratio_closehits_casc", ";ratio_closehits_casc", 400, 0,1);
  TH1D C_h_redToT_cascmuon("C_h_redToT_cascmuon", ";redToT_cascmuon", 400, 0, 4000);
  TH1D C_h_redToT_casc("C_h_redToT_casc", ";redToT_casc", 400, 0, 4000);
  TH1D C_h_myratio50_cascmuon("C_h_myratio50_cascmuon", ";myratio50_cascmuon", 400, 0, 1);
  TH1D C_h_myratio50_casc("C_h_myratio50_casc", ";myratio50_casc", 400, 0, 1);
  TH1D C_h_myratio30_cascmuon("C_h_myratio30_cascmuon", ";myratio30_cascmuon", 400, 0, 1);
  TH1D C_h_myratio30_muon("C_h_myratio30_muon", ";myratio30_muon", 400, 0, 1);
  TH1D C_h_myratio30_casc("C_h_myratio30_casc", ";myratio30_casc", 400, 0, 1);
  TH1D C_h_min_diff_sollik("C_h_min_diff_sollik",";min_diff_sollik", 400, 0, 4000);
  TH1D C_h_beta0_shower_deg("C_h_beta0_shower_deg", ";beta0_shower_deg", 400, 0, 4000);
  TH1D C_h_lik_shower("C_h_lik_shower", ";lik_shower", 400, 0, 4000);
  TH1D C_h_best_trk_pos_shower_z("C_h_best_trk_pos_shower_z", ";best_trk_pos_shower_z", 400, 0, 4000);
  TH1D C_h_normdlik("C_h_normdlik", ";normdlik", 400, 0, 4000);
  TH1D C_h_upsol("C_h_upsol", ";upsol", 100, 0, 60);
  TH1D C_h_Slen("C_h_Slen", ";Slen", 1000, 0, 1000);
  TH1D C_h_max_diff_sollik("C_h_max_diff_sollik", ";max_diff_sollik", 400, 0, 4000);
  TH1D C_h_diffangle_track_shower("C_h_diffangle_track_shower", ";diffangle_track_shower", 400, 0, 4000);
  TH1D C_h_myratio50_cascmuon_over_mu("C_h_myratio50_cascmuon_over_mu", ";myratio50_cascmuon_over_mu",400, 0, 1);
  TH1D C_h_myratio30_cascmuon_over_mu("C_h_myratio30_cascmuon_over_mu", ";myratio30_cascmuon_over_mu",400, 0, 1);
  TH1D C_h_myratio50_casc_over_mu("C_h_myratio50_casc_over_mu", ";myratio50_casc_over_mu",400, 0, 1);
  TH1D C_h_myratio30_casc_over_mu("C_h_myratio30_casc_over_mu", ";myratio30_casc_over_mu",400, 0, 1);
  TH1D C_h_ratio_closehits_cascmuon_over_mu("C_h_ratio_closehits_cascmuon_over_mu", ";ratio_closehits_cascmuon_over_mu",400, 0, 1);
  TH1D C_h_ratio_closehits_casc_over_mu("C_h_ratio_closehits_casc_over_mu", ";ratio_closehits_casc_over_mu",400, 0, 1);
  TH1D C_h_redToT_cascmuon_over_mu("C_h_redToT_cascmuon_over_mu", ";redToT_cascmuon_over_mu",400, 0, 1);
  TH1D C_h_redToT_casc_over_mu("C_h_redToT_casc_over_mu", ";redToT_casc_over_mu",400, 0, 1);
  TH1D C_h_bjorken_y("C_h_bjorken_y", ";bjorken_y",400, 0, 1);
  TH1D C_h_diffangle("C_h_diffangle", ";diffangle",180, 0, 180);
  TH1D C_h_LogE_mu ("C_h_LogE_mu","; LogE_mu ",100,0,10);
  TH1D C_h_LogEbundle ("C_h_LogEbundle","; LogEbundle ",100,0,10);
  //TH1Dc_//h_logE_mu_max ("//c_//h_logE_mu_max","; LogE_mu_max ",100,0,10);
  TH1D C_h_cos_zen_mu ("C_h_cos_zen_mu","; Cos(zenitC_h_mu)",100,-1,1);
  
  
  TH2D C_Lik_vs_NNhits ("C_Lik_vs_NNhits",";-Lik;NNhits",1000,0,1000,700,-100,600);
  TH2D C_Lik_vs_Q1value ("C_Lik_vs_Q1value","; ",500,-5,5,700,-100,600);
  TH2D C_NNhits_vs_Q1value ("C_NNhits_vs_Q1value","; ",500,-5,5,1000,0,1000);
  TH2D C_Nlines_vs_zenith ("C_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D C_Chere_Nlines_vs_zenith ("C_Chere_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D C_TrLengthIT_2_vs_Slen ("C_TrLengthIT_2_vs_Slen",";Slen;TrLengthIT_2",500,0,1000,500,0,1000);
  TH2D C_Diffangle_vs_run ("C_Diffangle_vs_run",";diffangle;run",4667,9333,14000,150,0,150);
  TH2D C_LogEreco_vs_CosZen ("C_LogEreco_vs_CosZen",";CosZen;LogEreco",100,-1,1,100,0,9);
  TH2D C_LogEreco_vs_LogEreco2 ("C_LogEreco_vs_LogEreco2",";LogEreco2;LogEreco",100,-1,10,100,-1,10);
  TH2D C_LogEreco_vs_LogE_mu ("C_LogEreco_vs_LogE_mu",";LogE_mu;LogEreco",100,-1,10,100,-1,10);
  TH2D C_LogEreco2_vs_LogE_mu ("C_LogEreco2_vs_LogE_mu",";LogE_mu;LogEreco2",100,-1,10,100,-1,10);
  TH2D C_Rate3Dmuon ("C_Rate3Dmuon","; run_nr",4667,9333,14000,10000,0,100);
  TH2D C_Rate3Dshower ("C_Rate3Dshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D C_RateMXshower ("C_RateMXshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D C_LogEreco_vs_LogE_neu ("C_LogEreco_vs_LogE_neu",";LogE_neu;LogEreco",100,-1,10,100,-1,10);
  TH2D C_LogEreco2_vs_LogE_neu ("C_LogEreco2_vs_LogE_neu",";LogE_neu;LogEreco2",100,-1,10,100,-1,10);
  TH2D C_LogEresolution_mu_vs_CosZen ("C_LogEresolution_mu_vs_CosZen",";CosZen;(LogEreco-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D C_LogEresolution_mu_cor_vs_CosZen ("C_LogEresolution_mu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D C_LogEresolution_neu_vs_CosZen ("C_LogEresolution_neu_vs_CosZen",";CosZen;(LogEreco-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  TH2D C_LogEresolution_neu_cor_vs_CosZen ("C_LogEresolution_neu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  
  //First cuts
  TH1D First_CosZen ("First_CosZen","First_CosZen; Cos(Zenith)",100,-1,1);
  TH1D First_Lik ("First_Lik","First_Lik; -Lik",700,-100,600);
  TH1D First_LogBeta0 ("First_LogBeta0","First_LogBeta0; Log(Beta0)",400,-4,0);
  TH1D First_Beta0_deg ("First_Beta0_deg","; LogBeta0deg",600,-4,2);
  TH1D First_DOMs ("First_DOMs","DOMs; # of DOMs",200,0,200);
  TH1D First_Cherenkov_DOMs ("First_Cherenkov_DOMs","Cherenkov_DOMs; # of Cherenkov DOMs",200,0,200);
  TH1D First_Q1value ("First_Q1value","First_Q1value; -Likelihood/Nhits",500,-5,5);
  TH1D First_NNhits ("First_NNhits","First_NNhits; # of Nhits",1000,0,1000);
  TH1D First_LogEreco ("First_LogEreco","; LogEreco ",100,0,8);
  TH1D First_TrLen  ("First_TrLen","; Track Length",1000,0,1000);
  TH1D First_Npe  ("First_Npe","; # of p.e",10000,0,10000);
  TH1D First_Ratiopmts  ("First_Ratiopmts","; Chere.PMTs/Reco.PMTs",100,0,1); 
  TH1D First_Ratiodoms  ("First_Ratiodoms","; Chere.DOMs/Reco.DOMs",100,0,1);
  TH2D First_Allreco_rate ("First_Allreco_rate","; run_nr",4667,9333,14000,10000,0,100);
  TH1D First_Rvalue ("First_Rvalue","First_R; R(m)",1000,0,1000);
  TH1D First_All_triggered_hits ("First_All_triggered_hits","First_All_triggered_hits; # of hits",1000,0,1000);
  TH1D First_Phi ("First_Phi","First_Phi; Azimuth [deg]",360,-180,180);
  TH1D First_Zenith_shower ("First_Zenith_shower","First_zenith_shower; [deg]",360,0,360);
  TH1D First_DLik ("First_DLik","; First_DLikelihood", 700,-100,600);
  TH1D First_Ratio6 ("First_Ratio6","; First_Ratio6",1200,0,1.2);
  TH1D First_MaxZenSol ("First_MaxZenSol","First_MaxZenSol; [deg]",360,0,360);
  TH1D First_RedToT_muon ("First_RedToT_muon","; Total_tot/nhits",300,0,300);
  TH1D First_Zenith ("First_Zenith","First_Zenith; [deg]",360,0,360);
  TH1D First_Delta_ZeniTh_track_shower ("First_Delta_ZeniTh_track_shower","First_Delta_ZeniTh_track_shower; [deg]",360,-180,180);
  TH1D First_MinZenSol ("First_MinZenSol","First_MinZenSol; [deg]",360,0,360);
  TH1D First_ItoLen ("First_ItoLen","; First_ItoLen",1200,0,1.2);
  TH1D First_Ratio330 ("First_Ratio330","; First_Ratio330",1200,0,1.2);
  TH1D First_Num_good_sol ("First_Num_good_sol",";# of good solutions",100,0,100);
  TH1D First_PosZ  ("First_PosZ","; First_PosZ",1000,0,1000);
  TH1D First_LogEnu ("First_LogEnu","; First_LogEnu ",100,0,8);
  TH1D First_LogAEnergy ("First_LogAEnergy","; First_LogAEnergy ",100,0,8);
  TH1D First_NLines ("First_NLines","; # of lines",12,0,12);
  TH1D First_LogEreco_cor ("First_LogEreco_cor","; LogEreco_cor ",100,0,8);
  TH1D First_Chere_NLines ("First_Chere_NLines","; # of lines",12,0,12);
  //TH1D First_Reco_Dclosest ("First_Reco_Dclosest","; d_closest ", 1000, 0, 1000);
  //TH1D First_Reco_tres ("First_Reco_tres","; reco_tres ", 1000, -500, 500);
  //TH1D First_Reco_Cos_Angle ("First_Reco_Cos_Angle","; reco_cos_angle ", 500, -360, 360);
  TH1D First_Ratio_upper_hemisphere ("First_Ratio_upper_hemisphere","; upper_hits/all_hits",100,0,1);
  TH1D First_Ratio_lower_hemisphere ("First_Ratio_lower_hemisphere","; lowerr_hits/all_hits",100,0,1);
  TH1D First_Ratio_border_upper_hemisphere ("First_Ratio_border_upper_hemisphere","; border_upper_hits/all_border_hits",100,0,1);
  TH1D First_Ratio_border_lower_hemisphere ("First_Ratio_border_lower_hemisphere","; border_lowerr_hits/all_border_hits",100,0,1);
  TH1D First_Ratio1 ("First_Ratio1","; Ratio1",1200,0,1.2);
  TH1D First_Ratio110 ("First_Ratio110","; Ratio110",1200,0,1.2);
  TH1D First_Ratio130 ("First_Ratio130","; Ratio130",1200,0,1.2);
  TH1D First_NtrackIT_ratio ("First_NtrackIT_ratio","; NtrackIT/Nhits",1200,0,1.2);
  TH1D First_TrLenIT_3  ("First_TrLenIT_3","; Track Length",1000,0,1000);

  TH1D  First_h_ToT_border_mu("First_h_ToT_border_mu",";ToT_border_mu ",400,0,4000);
  TH1D  First_h_ToT_border_casc("First_h_ToT_border_casc",";ToT_border_casc ",400,0,4000);
  TH1D  First_h_ToT_border_cascmu("First_h_ToT_border_cascmu",";ToT_border_cascmu ",400,0,4000);
  TH1D  First_h_ToT_trig("First_h_ToT_trig",";ToT_trig ",400,0,4000);
  TH1D  First_h_max_ToT_trig("First_h_max_ToT_trig",";max_ToT_trig ",400,0,4000);
  TH1D  First_h_ToT_IT("First_h_ToT_IT",";ToT_IT ",400,0,4000);
  TH1D  First_h_ToT_allIT("First_h_ToT_allIT",";ToT_allIT ",400,0,4000);
  TH1D  First_h_sum_ToT_casc("First_h_sum_ToT_casc",";sum_ToT_casc" ,400,0,4000);
  TH1D  First_h_sum_ToT_mu("First_h_sum_ToT_mu",";sum_ToT_mu" ,400,0,4000);
  TH1D  First_h_sum_ToT_casc_mu("First_h_sum_ToT_casc_mu",";sum_ToT_casc_mu" ,400,0,4000);
  //----------------------------------------------------------
  TH1D First_h_Ratio430 ("First_h_Ratio430","; Ratio430",1200,0,1.2);
  TH1D First_h_Ratio5 ("First_h_Ratio5","; Ratio5",1200,0,1.2);
  TH1D First_h_DiffDistance_mu  ("First_h_DiffDistance_mu","; dmax-dmin [m]",1400,0,1400);
  TH1D First_h_RatioCloseHits_mu  ("First_h_RatioCloseHits_mu","; nhits_400m/nhits",400,0,1);
  TH1D First_h_NtrackEarly("First_h_NtrackEarly", ";NtrackEarly",400,0,4000);
  TH1D First_h_D("First_h_D", ";D", 400, 0, 4000);
  TH1D First_h_NtrackIT30("First_h_NtrackIT30", ";NtrackIT30", 400, 0, 4000);
  TH1D First_h_NtrackLate("First_h_NtrackLate", ";NtrackLate", 400, 0, 4000);
  TH1D First_h_myratio50_muon("First_h_myratio50_muon", ";myratio50_muon", 400, 0, 1);
  TH1D First_h_ratio_cherenkov_lines("First_h_ratio_cherenkov_lines", ";ratio_cherenkov_lines", 400, 0, 1);
  TH1D First_h_SnpeT("First_h_SnpeT", ";SnpeT", 400, 0, 4000);
  TH1D First_h_mean_tres_it("First_h_mean_tres_it", ";mean_tres_it", 400, 0, 4000);
  TH1D First_h_max_lik_up("First_h_max_lik_up", ";max_lik_up", 400, 0, 4000);
  TH1D First_h_max_lik_down("First_h_max_lik_down", ";max_lik_down", 400, 0, 4000);
  TH1D First_h_diff_theta("First_h_diff_theta", ";diff_theta", 400, 0, 4000);
  TH1D First_h_diff_dist_casc_mu("First_h_diff_dist_casc_mu", ";diff_dist_casc_mu", 400, 0, 4000);
  TH1D First_h_diff_dist_casc("First_h_diff_dist_casc", ";diff_dist_casc", 400, 0, 4000);
  TH1D First_h_ratio_closehits_cascmuon("First_h_ratio_closehits_cascmuon", ";ratio_closehits_cascmuon", 400, 0, 1);
  TH1D First_h_ratio_closehits_casc("First_h_ratio_closehits_casc", ";ratio_closehits_casc", 400, 0,1);
  TH1D First_h_redToT_cascmuon("First_h_redToT_cascmuon", ";redToT_cascmuon", 400, 0, 4000);
  TH1D First_h_redToT_casc("First_h_redToT_casc", ";redToT_casc", 400, 0, 4000);
  TH1D First_h_myratio50_cascmuon("First_h_myratio50_cascmuon", ";myratio50_cascmuon", 400, 0, 1);
  TH1D First_h_myratio50_casc("First_h_myratio50_casc", ";myratio50_casc", 400, 0, 1);
  TH1D First_h_myratio30_cascmuon("First_h_myratio30_cascmuon", ";myratio30_cascmuon", 400, 0, 1);
  TH1D First_h_myratio30_muon("First_h_myratio30_muon", ";myratio30_muon", 400, 0, 1);
  TH1D First_h_myratio30_casc("First_h_myratio30_casc", ";myratio30_casc", 400, 0, 1);
  TH1D First_h_min_diff_sollik("First_h_min_diff_sollik",";min_diff_sollik", 400, 0, 4000);
  TH1D First_h_beta0_shower_deg("First_h_beta0_shower_deg", ";beta0_shower_deg", 400, 0, 4000);
  TH1D First_h_lik_shower("First_h_lik_shower", ";lik_shower", 400, 0, 4000);
  TH1D First_h_best_trk_pos_shower_z("First_h_best_trk_pos_shower_z", ";best_trk_pos_shower_z", 400, 0, 4000);
  TH1D First_h_normdlik("First_h_normdlik", ";normdlik", 400, 0, 4000);
  TH1D First_h_upsol("First_h_upsol", ";upsol", 100, 0, 60);
  TH1D First_h_Slen("First_h_Slen", ";Slen", 1000, 0, 1000);
  TH1D First_h_max_diff_sollik("First_h_max_diff_sollik", ";max_diff_sollik", 400, 0, 4000);
  TH1D First_h_diffangle_track_shower("First_h_diffangle_track_shower", ";diffangle_track_shower", 400, 0, 4000);
  TH1D First_h_myratio50_cascmuon_over_mu("First_h_myratio50_cascmuon_over_mu", ";myratio50_cascmuon_over_mu",400, 0, 1);
  TH1D First_h_myratio30_cascmuon_over_mu("First_h_myratio30_cascmuon_over_mu", ";myratio30_cascmuon_over_mu",400, 0, 1);
  TH1D First_h_myratio50_casc_over_mu("First_h_myratio50_casc_over_mu", ";myratio50_casc_over_mu",400, 0, 1);
  TH1D First_h_myratio30_casc_over_mu("First_h_myratio30_casc_over_mu", ";myratio30_casc_over_mu",400, 0, 1);
  TH1D First_h_ratio_closehits_cascmuon_over_mu("First_h_ratio_closehits_cascmuon_over_mu", ";ratio_closehits_cascmuon_over_mu",400, 0, 1);
  TH1D First_h_ratio_closehits_casc_over_mu("First_h_ratio_closehits_casc_over_mu", ";ratio_closehits_casc_over_mu",400, 0, 1);
  TH1D First_h_redToT_cascmuon_over_mu("First_h_redToT_cascmuon_over_mu", ";redToT_cascmuon_over_mu",400, 0, 1);
  TH1D First_h_redToT_casc_over_mu("First_h_redToT_casc_over_mu", ";redToT_casc_over_mu",400, 0, 1);
  TH1D First_h_bjorken_y("First_h_bjorken_y", ";bjorken_y",400, 0, 1);
  TH1D First_h_diffangle("First_h_diffangle", ";diffangle",180, 0, 180);
  TH1D First_h_LogE_mu ("First_h_LogE_mu","; LogE_mu ",100,0,10);
  TH1D First_h_LogEbundle ("First_h_LogEbundle","; LogEbundle ",100,0,10);
  //TH1DFirst_//h_logE_mu_max ("//First_//h_logE_mu_max","; LogE_mu_max ",100,0,10);
  TH1D First_h_cos_zen_mu ("First_h_cos_zen_mu","; Cos(zenitFirst_h_mu)",100,-1,1);

  TH2D First_Lik_vs_NNhits ("First_Lik_vs_NNhits","; ",1000,0,1000,700,-100,600);
  TH2D First_Lik_vs_Q1value ("First_Lik_vs_Q1value","; ",500,-5,5,700,-100,600);
  TH2D First_NNhits_vs_Q1value ("First_NNhits_vs_Q1value","; ",500,-5,5,1000,0,1000);
  TH2D First_Nlines_vs_zenith ("First_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D First_Chere_Nlines_vs_zenith ("First_Chere_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D First_TrLengthIT_2_vs_Slen ("First_TrLengthIT_2_vs_Slen",";Slen;TrLengthIT_2",500,0,1000,500,0,1000);
  TH2D First_Diffangle_vs_run ("First_Diffangle_vs_run",";diffangle;run",4667,9333,14000,150,0,150);
  TH2D First_LogEreco_vs_CosZen ("First_LogEreco_vs_CosZen",";CosZen;LogEreco",100,-1,1,100,0,9);
  TH2D First_LogEreco_vs_LogEreco2 ("First_LogEreco_vs_LogEreco2",";LogEreco2;LogEreco",100,-1,10,100,-1,10);
  TH2D First_LogEreco_vs_LogE_mu ("First_LogEreco_vs_LogE_mu",";LogE_mu;LogEreco",100,-1,10,100,-1,10);
  TH2D First_LogEreco2_vs_LogE_mu ("First_LogEreco2_vs_LogE_mu",";LogE_mu;LogEreco2",100,-1,10,100,-1,10);
  TH2D First_Rate3Dmuon ("First_Rate3Dmuon","; run_nr",4667,9333,14000,10000,0,100);
  TH2D First_Rate3Dshower ("First_Rate3Dshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D First_RateMXshower ("First_RateMXshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D First_LogEreco_vs_LogE_neu ("First_LogEreco_vs_LogE_neu",";LogE_neu;LogEreco",100,-1,10,100,-1,10);
  TH2D First_LogEreco2_vs_LogE_neu ("First_LogEreco2_vs_LogE_neu",";LogE_neu;LogEreco2",100,-1,10,100,-1,10);
  TH2D First_LogEresolution_mu_vs_CosZen ("First_LogEresolution_mu_vs_CosZen",";CosZen;(LogEreco-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D First_LogEresolution_mu_cor_vs_CosZen ("First_LogEresolution_mu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D First_LogEresolution_neu_vs_CosZen ("First_LogEresolution_neu_vs_CosZen",";CosZen;(LogEreco-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  TH2D First_LogEresolution_neu_cor_vs_CosZen ("First_LogEresolution_neu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  
  //-------------Same but for up going
  TH1D Up_CosZen ("Up_CosZen","; Cos(reco_zenith)",100,-1,1);
  TH1D Up_Lik ("Up_Lik","; Likelihood", 700,-100,600);
  TH1D Up_LogBeta0 ("Up_LogBeta0","; LogBeta0",400,-4,0);
  TH1D Up_Beta0_deg ("Up_Beta0_deg","; LogBeta0deg",600,-4,2);
  TH1D Up_DOMs ("Up_DOMs","DOMs; # of DOMs",200,0,200);
  TH1D Up_Cherenkov_DOMs ("Up_Cherenkov_DOMs","Cherenkov_DOMs; # of Cherenkov DOMs",200,0,200);
  TH1D Up_Q1value ("Up_Q1value","Up_Q1value; -Likelihood/Nhits",500,-5,5);
  TH1D Up_NNhits ("Up_NNhits","Up_NNhits; # of Nhits",1000,0,1000); 
  TH1D Up_LogEreco ("Up_LogEreco","; LogEreco ",100,0,8);
  TH1D Up_TrLen  ("Up_TrLen","; Track Length",1000,0,1000);
  TH1D Up_Npe  ("Up_Npe","; # of p.e",10000,0,10000);
  TH1D Up_Ratiopmts  ("Up_Ratiopmts","; Chere.PMTs/Reco.PMTs",100,0,1); 
  TH1D Up_Ratiodoms  ("Up_Ratiodoms","; Chere.DOMs/Reco.DOMs",100,0,1);
  TH2D Up_Allreco_rate ("Up_Allreco_rate","; run_nr",4667,9333,14000,10000,0,100);
  TH1D Up_Rvalue ("Up_Rvalue","Up_R; R(m)",1000,0,1000);
  TH1D Up_All_triggered_hits ("Up_All_triggered_hits","Up_All_triggered_hits; # of hits",1000,0,1000);
  TH1D Up_Phi ("Up_Phi","Up_Phi; Azimuth [deg]",360,-180,180);
  TH1D Up_Zenith_shower ("Up_Zenith_shower","Up_zenith_shower; [deg]",360,0,360);
  TH1D Up_DLik ("Up_DLik","; Up_DLikelihood", 700,-100,600);
  TH1D Up_Ratio6 ("Up_Ratio6","; Up_Ratio6",1200,0,1.2);
  TH1D Up_MaxZenSol ("Up_MaxZenSol","Up_MaxZenSol; [deg]",360,0,360);
  TH1D Up_RedToT_muon ("Up_RedToT_muon","; Total_tot/nhits",300,0,300);
  TH1D Up_Zenith ("Up_Zenith","Up_Zenith; [deg]",360,0,360);
  TH1D Up_Delta_ZeniTh_track_shower ("Up_Delta_ZeniTh_track_shower","Up_Delta_ZeniTh_track_shower; [deg]",360,-180,180);
  TH1D Up_MinZenSol ("Up_MinZenSol","Up_MinZenSol; [deg]",360,0,360);
  TH1D Up_ItoLen ("Up_ItoLen","; Up_ItoLen",1200,0,1.2);
  TH1D Up_Ratio330 ("Up_Ratio330","; Up_Ratio330",1200,0,1.2);
  TH1D Up_Num_good_sol ("Up_Num_good_sol",";# of good solutions",100,0,100);
  TH1D Up_PosZ  ("Up_PosZ","; Up_PosZ",1000,0,1000);
  TH1D Up_LogEnu ("Up_LogEnu","; Up_LogEnu ",100,0,8);
  TH1D Up_LogAEnergy ("Up_LogAEnergy","; Up_LogAEnergy ",100,0,8);
  TH1D Up_NLines ("Up_NLines","; # of lines",12,0,12);
  TH1D Up_LogEreco_cor ("Up_LogEreco_cor","; LogEreco_cor ",100,0,8);
  TH1D Up_Chere_NLines ("Up_Chere_NLines","; # of lines",12,0,12);
  //TH1D Up_Reco_Dclosest ("Up_Reco_Dclosest","; d_closest ", 1000, 0, 1000);
  //TH1D Up_Reco_tres ("Up_Reco_tres","; reco_tres ", 1000, -500, 500);
  //TH1D Up_Reco_Cos_Angle ("Up_Reco_Cos_Angle","; reco_cos_angle ", 500, -360, 360);
  TH1D Up_Ratio_upper_hemisphere ("Up_Ratio_upper_hemisphere","; upper_hits/all_hits",100,0,1);
  TH1D Up_Ratio_lower_hemisphere ("Up_Ratio_lower_hemisphere","; lowerr_hits/all_hits",100,0,1);
  TH1D Up_Ratio_border_upper_hemisphere ("Up_Ratio_border_upper_hemisphere","; border_upper_hits/all_border_hits",100,0,1);
  TH1D Up_Ratio_border_lower_hemisphere ("Up_Ratio_border_lower_hemisphere","; border_lowerr_hits/all_border_hits",100,0,1);
  TH1D Up_Ratio1 ("Up_Ratio1","; Ratio1",1200,0,1.2);
  TH1D Up_Ratio110 ("Up_Ratio110","; Ratio110",1200,0,1.2);
  TH1D Up_Ratio130 ("Up_Ratio130","; Ratio130",1200,0,1.2);
  TH1D Up_NtrackIT_ratio ("Up_NtrackIT_ratio","; NtrackIT/Nhits",1200,0,1.2);
  TH1D Up_TrLenIT_3  ("Up_TrLenIT_3","; Track Length",1000,0,1000);

  TH1D Up_h_ToT_border_mu("Up_h_ToT_border_mu",";Up_ToT_border_mu ",400,0,4000);
  TH1D Up_h_ToT_border_casc("Up_h_ToT_border_casc",";Up_ToT_border_casc ",400,0,4000);
  TH1D Up_h_ToT_border_cascmu("Up_h_ToT_border_cascmu",";Up_ToT_border_cascmu ",400,0,4000);
  TH1D Up_h_ToT_trig("Up_h_ToT_trig",";Up_ToT_trig ",400,0,4000);
  TH1D Up_h_max_ToT_trig("h_max_Up_ToT_trig",";max_Up_ToT_trig ",400,0,4000);
  TH1D Up_h_ToT_IT("Up_h_ToT_IT",";Up_ToT_IT ",400,0,4000);
  TH1D Up_h_ToT_allIT("Up_h_ToT_allIT",";Up_ToT_allIT ",400,0,4000);
  TH1D Up_h_sum_ToT_casc("Up_h_sum_ToT_casc",";Up_sum_ToT_casc" ,400,0,4000);
  TH1D Up_h_sum_ToT_mu("Up_h_sum_ToT_mu",";Up_sum_ToT_mu" ,400,0,4000);
  TH1D Up_h_sum_ToT_casc_mu("Up_h_sum_ToT_casc_mu",";Up_sum_ToT_casc_mu" ,400,0,4000);

  //---------------------------------------------------
  TH1D Up_h_Ratio430 ("Up_h_Ratio430","; Ratio430",1200,0,1.2);
  TH1D Up_h_Ratio5 ("Up_h_Ratio5","; Ratio5",1200,0,1.2);
  TH1D Up_h_DiffDistance_mu  ("Up_h_DiffDistance_mu","; dmax-dmin [m]",12000,0,12000);
  TH1D Up_h_RatioCloseHits_mu  ("Up_h_RatioCloseHits_mu","; nhits_400m/nhits",400,0,1);
  TH1D Up_h_NtrackEarly("Up_h_NtrackEarly", ";NtrackEarly",400,0,4000);
  TH1D Up_h_D("Up_h_D"," ;D", 400, 0, 4000);
  TH1D Up_h_NtrackIT30("Up_h_NtrackIT30", ";NtrackIT30", 400, 0, 4000);
  TH1D Up_h_NtrackLate("Up_h_NtrackLate", ";Up_NtrackLate", 400, 0, 4000);
  TH1D Up_h_myratio50_muon("Up_h_myratio50_muon", ";Up_myratio50_muon", 400, 0, 4000);
  TH1D Up_h_ratio_cherenkov_lines("Up_h_ratio_cherenkov_lines", ";Up_ratio_cherenkov_lines", 400, 0, 4000);
  TH1D Up_h_SnpeT("Up_h_SnpeT", ";Up_SnpeT", 400, 0, 4000);
  TH1D Up_h_mean_tres_it("Up_h_mean_tres_it", ";Up_mean_tres_it", 400, 0, 4000);
  TH1D Up_h_max_lik_up("Up_h_max_lik_up", ";Up_max_lik_up", 400, 0, 4000);
  TH1D Up_h_max_lik_down("Up_h_max_lik_down", ";Up_max_lik_down", 400, 0, 4000);
  TH1D Up_h_diff_theta("Up_h_diff_theta", ";Up_diff_theta", 400, 0, 4000);
  TH1D Up_h_diff_dist_casc_mu("Up_h_diff_dist_casc_mu", ";Up_diff_dist_casc_mu", 400, 0, 4000);
  TH1D Up_h_diff_dist_casc("Up_h_diff_dist_casc", ";Up_diff_dist_casc", 400, 0, 4000);
  TH1D Up_h_ratio_closehits_cascmuon("Up_h_ratio_closehits_cascmuon", ";Up_ratio_closehits_cascmuon", 400, 0, 4000);
  TH1D Up_h_ratio_closehits_casc("Up_h_ratio_closehits_casc", ";Up_ratio_closehits_casc", 400, 0, 4000);
  TH1D Up_h_redToT_cascmuon("Up_h_redToT_cascmuon", ";Up_redToT_cascmuon", 400, 0, 4000);
  TH1D Up_h_redToT_casc("Up_h_redToT_casc", ";Up_redToT_casc", 400, 0, 4000);
  TH1D Up_h_myratio50_cascmuon("Up_h_myratio50_cascmuon", ";Up_myratio50_cascmuon", 400, 0, 4000);
  TH1D Up_h_myratio50_casc("Up_h_myratio50_casc", ";Up_myratio50_casc", 400, 0, 4000);
  TH1D Up_h_myratio30_cascmuon("Up_h_myratio30_cascmuon", ";Up_myratio30_cascmuon", 400, 0, 4000);
  TH1D Up_h_myratio30_muon("Up_h_myratio30_muon", ";Up_myratio30_muon", 400, 0, 4000);
  TH1D Up_h_myratio30_casc("Up_h_myratio30_casc", ";Up_myratio30_casc", 400, 0, 4000);
  TH1D Up_h_min_diff_sollik("Up_h_min_diff_sollik",";Up_min_diff_sollik", 400, 0, 4000);
  TH1D Up_h_beta0_shower_deg("Up_h_beta0_shower_deg", ";Up_beta0_shower_deg", 400, 0, 4000);
  TH1D Up_h_lik_shower("Up_h_lik_shower", ";Up_lik_shower", 400, 0, 4000);
  TH1D Up_h_best_trk_pos_shower_z("Up_h_best_trk_pos_shower_z", ";Up_best_trk_pos_shower_z", 400, 0, 4000);
  TH1D Up_h_normdlik("Up_h_normdlik", ";Up_normdlik", 400, 0, 4000);
  TH1D Up_h_upsol("Up_h_upsol" ,";Up_upsol", 100, 0, 60);
  TH1D Up_h_Slen("Up_h_Slen", ";Up_Slen", 1000, 0, 1000);
  TH1D Up_h_max_diff_sollik("Up_h_max_diff_sollik", ";Up_max_diff_sollik", 400, 0, 4000);
  TH1D Up_h_diffangle_track_shower("Up_h_diffangle_track_shower", ";Up_diffangle_track_shower", 400, 0, 4000);
  TH1D Up_h_myratio50_cascmuon_over_mu("Up_h_myratio50_cascmuon_over_mu", ";myratio50_cascmuon_over_mu",400, 0, 1);
  TH1D Up_h_myratio30_cascmuon_over_mu("Up_h_myratio30_cascmuon_over_mu", ";myratio30_cascmuon_over_mu",400, 0, 1);
  TH1D Up_h_myratio50_casc_over_mu("Up_h_myratio50_casc_over_mu", ";myratio50_casc_over_mu",400, 0, 1);
  TH1D Up_h_myratio30_casc_over_mu("Up_h_myratio30_casc_over_mu", ";myratio30_casc_over_mu",400, 0, 1);
  TH1D Up_h_ratio_closehits_cascmuon_over_mu("Up_h_ratio_closehits_cascmuon_over_mu", ";ratio_closehits_cascmuon_over_mu",400, 0, 1);
  TH1D Up_h_ratio_closehits_casc_over_mu("Up_h_ratio_closehits_casc_over_mu", ";ratio_closehits_casc_over_mu",400, 0, 1);
  TH1D Up_h_redToT_cascmuon_over_mu("Up_h_redToT_cascmuon_over_mu", ";redToT_cascmuon_over_mu",400, 0, 1);
  TH1D Up_h_redToT_casc_over_mu("Up_h_redToT_casc_over_mu", ";redToT_casc_over_mu",400, 0, 1);
  TH1D Up_h_bjorken_y("Up_h_bjorken_y", ";bjorken_y",400, 0, 1);
  TH1D Up_h_diffangle("Up_h_diffangle", ";diffangle",180, 0, 180);
  TH1D Up_h_LogE_mu ("Up_h_LogE_mu","; LogE_mu ",100,0,10);
  TH1D Up_h_LogEbundle ("Up_h_LogEbundle","; LogEbundle ",100,0,10);
  //TH1D Up_//h_logE_mu_max ("Up_//h_logE_mu_max","; LogE_mu_max ",100,0,10);
  TH1D Up_h_cos_zen_mu ("Up_h_cos_zen_mu","; Cos(zenitUp_h_mu)",100,-1,1);
  //----------------------------------------------------------------

  TH2D Up_Lik_vs_NNhits ("Up_Lik_vs_NNhits","; ",1000,0,1000,700,-100,600);
  TH2D Up_Lik_vs_Q1value ("Up_Lik_vs_Q1value","; ",500,-5,5,700,-100,600);
  TH2D Up_NNhits_vs_Q1value ("Up_NNhits_vs_Q1value","; ",500,-5,5,1000,0,1000);
  TH2D Up_Nlines_vs_zenith ("Up_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D Up_Chere_Nlines_vs_zenith ("Up_Chere_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D Up_TrLengthIT_2_vs_Slen ("Up_TrLengthIT_2_vs_Slen",";Slen;TrLengthIT_2",500,0,1000,500,0,1000);
  TH2D Up_Diffangle_vs_run ("Up_Diffangle_vs_run",";diffangle;run",4667,9333,14000,150,0,150);
  TH2D Up_LogEreco_vs_CosZen ("Up_LogEreco_vs_CosZen",";CosZen;LogEreco",100,-1,1,100,0,9);
  TH2D Up_LogEreco_vs_LogEreco2 ("Up_LogEreco_vs_LogEreco2",";LogEreco2;LogEreco",100,-1,10,100,-1,10);
  TH2D Up_LogEreco_vs_LogE_mu ("Up_LogEreco_vs_LogE_mu",";LogE_mu;LogEreco",100,-1,10,100,-1,10);
  TH2D Up_LogEreco2_vs_LogE_mu ("Up_LogEreco2_vs_LogE_mu",";LogE_mu;LogEreco2",100,-1,10,100,-1,10);
  TH2D Up_Rate3Dmuon ("Up_Rate3Dmuon","; run_nr",4667,9333,14000,10000,0,100);
  TH2D Up_Rate3Dshower ("Up_Rate3Dshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D Up_RateMXshower ("Up_RateMXshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D Up_LogEreco_vs_LogE_neu ("Up_LogEreco_vs_LogE_neu",";LogE_neu;LogEreco",100,-1,10,100,-1,10);
  TH2D Up_LogEreco2_vs_LogE_neu ("Up_LogEreco2_vs_LogE_neu",";LogE_neu;LogEreco2",100,-1,10,100,-1,10);
  TH2D Up_LogEresolution_mu_vs_CosZen ("Up_LogEresolution_mu_vs_CosZen",";CosZen;(LogEreco-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D Up_LogEresolution_mu_cor_vs_CosZen ("Up_LogEresolution_mu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D Up_LogEresolution_neu_vs_CosZen ("Up_LogEresolution_neu_vs_CosZen",";CosZen;(LogEreco-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  TH2D Up_LogEresolution_neu_cor_vs_CosZen ("Up_LogEresolution_neu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);

  //Anti-noise cuts
  TH1D C_Up_CosZen ("C_Up_CosZen","; Cos(reco_zenith)",100,-1,1);
  TH1D C_Up_Lik ("C_Up_Lik","; Likelihood", 700,-100,600);
  TH1D C_Up_LogBeta0 ("C_Up_LogBeta0","; LogBeta0",400,-4,0);
  TH1D C_Up_Beta0_deg ("C_Up_Beta0_deg","; LogBeta0deg",600,-4,2);
  TH1D C_Up_DOMs ("C_Up_DOMs","C_Up_DOMs; # of DOMs",200,0,200);
  TH1D C_Up_Cherenkov_DOMs ("C_Up_Cherenkov_DOMs","C_Up_Cherenkov_DOMs; # of Cherenkov DOMs",200,0,200);
  TH1D C_Up_Q1value ("C_Up_Q1value","C_Up_Q1value; -Likelihood/Nhits",500,-5,5);
  TH1D C_Up_NNhits ("C_Up_NNhits","C_Up_NNhits; # of Nhits",1000,0,1000); 
  TH1D C_Up_LogEreco ("C_Up_LogEreco","; LogEreco ",100,0,8);
  TH1D C_Up_TrLen  ("C_Up_TrLen","; Track Length",1000,0,1000);
  TH1D C_Up_Npe  ("C_Up_Npe","; # of p.e",10000,0,10000);
  TH1D C_Up_Ratiopmts  ("C_Up_Ratiopmts","; Chere.PMTs/Reco.PMTs",100,0,1); 
  TH1D C_Up_Ratiodoms  ("C_Up_Ratiodoms","; Chere.DOMs/Reco.DOMs",100,0,1);
  TH2D C_Up_Allreco_rate ("C_Up_Allreco_rate","; run_nr",4667,9333,14000,10000,0,100);
  TH1D C_Up_Rvalue ("C_Up_Rvalue","Up_R; R(m)",1000,0,1000);
  TH1D C_Up_All_triggered_hits ("C_Up_All_triggered_hits","C_Up_All_triggered_hits; # of hits",1000,0,1000);
  TH1D C_Up_Phi ("C_Up_Phi","C_Up_Phi; Azimuth [deg]",360,-180,180);
  TH1D C_Up_Zenith_shower ("C_Up_Zenith_shower","C_Up_zenith_shower; [deg]",360,0,360);
  TH1D C_Up_DLik ("C_Up_DLik","; C_Up_DLikelihood", 700,-100,600);
  TH1D C_Up_Ratio6 ("C_Up_Ratio6","; C_Up_Ratio6",1200,0,1.2);
  TH1D C_Up_MaxZenSol ("C_Up_MaxZenSol","C_Up_MaxZenSol; [deg]",360,0,360);
  TH1D C_Up_RedToT_muon ("C_Up_RedToT_muon","; Total_tot/nhits",300,0,300);
  TH1D C_Up_Zenith ("C_Up_Zenith","C_Up_Zenith; [deg]",360,0,360);
  TH1D C_Up_Delta_ZeniTh_track_shower ("C_Up_Delta_ZeniTh_track_shower","C_Up_Delta_ZeniTh_track_shower; [deg]",360,-180,180);
  TH1D C_Up_MinZenSol ("C_Up_MinZenSol","C_Up_MinZenSol; [deg]",360,0,360);
  TH1D C_Up_ItoLen ("C_Up_ItoLen","; C_Up_ItoLen",1200,0,1.2);
  TH1D C_Up_Ratio330 ("C_Up_Ratio330","; C_Up_Ratio330",1200,0,1.2);
  TH1D C_Up_Num_good_sol ("C_Up_Num_good_sol",";# of good solutions",100,0,100);
  TH1D C_Up_PosZ  ("C_Up_PosZ","; C_Up_PosZ",1000,0,1000);
  TH1D C_Up_LogEnu ("C_Up_LogEnu","; C_Up_LogEnu ",100,0,8);
  TH1D C_Up_LogAEnergy ("C_Up_LogAEnergy","; C_Up_LogAEnergy ",100,0,8);
  TH1D C_Up_NLines ("C_Up_NLines","; # of lines",12,0,12);
  TH1D C_Up_LogEreco_cor ("C_Up_LogEreco_cor","; LogEreco_cor ",100,0,8);
  TH1D C_Up_Chere_NLines ("C_Up_Chere_NLines","; # of lines",12,0,12);
  //TH1D C_Up_Reco_Dclosest ("C_Up_Reco_Dclosest","; d_closest ", 1000, 0, 1000);
  //TH1D C_Up_Reco_tres ("C_Up_Reco_tres","; reco_tres ", 1000, -500, 500);
  //TH1D C_Up_Reco_Cos_Angle ("C_Up_Reco_Cos_Angle","; reco_cos_angle ", 500, -360, 360);
  TH1D C_Up_Ratio_upper_hemisphere ("C_Up_Ratio_upper_hemisphere","; upper_hits/all_hits",100,0,1);
  TH1D C_Up_Ratio_lower_hemisphere ("C_Up_Ratio_lower_hemisphere","; lowerr_hits/all_hits",100,0,1);
  TH1D C_Up_Ratio_border_upper_hemisphere ("C_Up_Ratio_border_upper_hemisphere","; border_upper_hits/all_border_hits",100,0,1);
  TH1D C_Up_Ratio_border_lower_hemisphere ("C_Up_Ratio_border_lower_hemisphere","; border_lowerr_hits/all_border_hits",100,0,1);
  TH1D C_Up_Ratio1 ("C_Up_Ratio1","; Ratio1",1200,0,1.2);
  TH1D C_Up_Ratio110 ("C_Up_Ratio110","; Ratio110",1200,0,1.2);
  TH1D C_Up_Ratio130 ("C_Up_Ratio130","; Ratio130",1200,0,1.2);
  TH1D C_Up_NtrackIT_ratio ("C_Up_NtrackIT_ratio","; NtrackIT/Nhits",1200,0,1.2);
  TH1D C_Up_TrLenIT_3  ("C_Up_TrLenIT_3","; Track Length",1000,0,1000);

  TH1D  C_Up_h_ToT_border_mu("C_Up_h_ToT_border_mu",";ToT_border_mu ",400,0,4000);
  TH1D  C_Up_h_ToT_border_casc("C_Up_h_ToT_border_casc",";ToT_border_casc ",400,0,4000);
  TH1D  C_Up_h_ToT_border_cascmu("C_Up_h_ToT_border_cascmu",";ToT_border_cascmu ",400,0,4000);
  TH1D  C_Up_h_ToT_trig("C_Up_h_ToT_trig",";ToT_trig ",400,0,4000);
  TH1D  C_Up_h_max_ToT_trig("C_Up_h_max_ToT_trig",";max_ToT_trig ",400,0,4000);
  TH1D  C_Up_h_ToT_IT("C_Up_h_ToT_IT",";ToT_IT ",400,0,4000);
  TH1D  C_Up_h_ToT_allIT("C_Up_h_ToT_allIT",";ToT_allIT ",400,0,4000);
  TH1D  C_Up_h_sum_ToT_casc("C_Up_h_sum_ToT_casc",";sum_ToT_casc" ,400,0,4000);
  TH1D  C_Up_h_sum_ToT_mu("C_Up_h_sum_ToT_mu",";sum_ToT_mu" ,400,0,4000);
  TH1D  C_Up_h_sum_ToT_casc_mu("C_Up_h_sum_ToT_casc_mu",";sum_ToT_casc_mu" ,400,0,4000);
  //----------------------------------------------------------
  TH1D C_Up_h_Ratio430 ("C_Up_h_Ratio430","; Ratio430",1200,0,1.2);
  TH1D C_Up_h_Ratio5 ("C_Up_h_Ratio5","; Ratio5",1200,0,1.2);
  TH1D C_Up_h_DiffDistance_mu  ("C_Up_h_DiffDistance_mu","; dmax-dmin [m]",1400,0,1400);
  TH1D C_Up_h_RatioCloseHits_mu  ("C_Up_h_RatioCloseHits_mu","; nhits_400m/nhits",400,0,1);
  TH1D C_Up_h_NtrackEarly("C_Up_h_NtrackEarly", ";NtrackEarly",400,0,4000);
  TH1D C_Up_h_D("C_Up_h_D", ";D", 400, 0, 4000);
  TH1D C_Up_h_NtrackIT30("C_Up_h_NtrackIT30", ";NtrackIT30", 400, 0, 4000);
  TH1D C_Up_h_NtrackLate("C_Up_h_NtrackLate", ";NtrackLate", 400, 0, 4000);
  TH1D C_Up_h_myratio50_muon("C_Up_h_myratio50_muon", ";myratio50_muon", 400, 0, 1);
  TH1D C_Up_h_ratio_cherenkov_lines("C_Up_h_ratio_cherenkov_lines", ";ratio_cherenkov_lines", 400, 0, 1);
  TH1D C_Up_h_SnpeT("C_Up_h_SnpeT", ";SnpeT", 400, 0, 4000);
  TH1D C_Up_h_mean_tres_it("C_Up_h_mean_tres_it", ";mean_tres_it", 400, 0, 4000);
  TH1D C_Up_h_max_lik_up("C_Up_h_max_lik_up", ";max_lik_up", 400, 0, 4000);
  TH1D C_Up_h_max_lik_down("C_Up_h_max_lik_down", ";max_lik_down", 400, 0, 4000);
  TH1D C_Up_h_diff_theta("C_Up_h_diff_theta", ";diff_theta", 400, 0, 4000);
  TH1D C_Up_h_diff_dist_casc_mu("C_Up_h_diff_dist_casc_mu", ";diff_dist_casc_mu", 400, 0, 4000);
  TH1D C_Up_h_diff_dist_casc("C_Up_h_diff_dist_casc", ";diff_dist_casc", 400, 0, 4000);
  TH1D C_Up_h_ratio_closehits_cascmuon("C_Up_h_ratio_closehits_cascmuon", ";ratio_closehits_cascmuon", 400, 0, 1);
  TH1D C_Up_h_ratio_closehits_casc("C_Up_h_ratio_closehits_casc", ";ratio_closehits_casc", 400, 0,1);
  TH1D C_Up_h_redToT_cascmuon("C_Up_h_redToT_cascmuon", ";redToT_cascmuon", 400, 0, 4000);
  TH1D C_Up_h_redToT_casc("C_Up_h_redToT_casc", ";redToT_casc", 400, 0, 4000);
  TH1D C_Up_h_myratio50_cascmuon("C_Up_h_myratio50_cascmuon", ";myratio50_cascmuon", 400, 0, 1);
  TH1D C_Up_h_myratio50_casc("C_Up_h_myratio50_casc", ";myratio50_casc", 400, 0, 1);
  TH1D C_Up_h_myratio30_cascmuon("C_Up_h_myratio30_cascmuon", ";myratio30_cascmuon", 400, 0, 1);
  TH1D C_Up_h_myratio30_muon("C_Up_h_myratio30_muon", ";myratio30_muon", 400, 0, 1);
  TH1D C_Up_h_myratio30_casc("C_Up_h_myratio30_casc", ";myratio30_casc", 400, 0, 1);
  TH1D C_Up_h_min_diff_sollik("C_Up_h_min_diff_sollik",";min_diff_sollik", 400, 0, 4000);
  TH1D C_Up_h_beta0_shower_deg("C_Up_h_beta0_shower_deg", ";beta0_shower_deg", 400, 0, 4000);
  TH1D C_Up_h_lik_shower("C_Up_h_lik_shower", ";lik_shower", 400, 0, 4000);
  TH1D C_Up_h_best_trk_pos_shower_z("C_Up_h_best_trk_pos_shower_z", ";best_trk_pos_shower_z", 400, 0, 4000);
  TH1D C_Up_h_normdlik("C_Up_h_normdlik", ";normdlik", 400, 0, 4000);
  TH1D C_Up_h_upsol("C_Up_h_upsol", ";upsol", 100, 0, 60);
  TH1D C_Up_h_Slen("C_Up_h_Slen", ";Slen", 1000, 0, 1000);
  TH1D C_Up_h_max_diff_sollik("C_Up_h_max_diff_sollik", ";max_diff_sollik", 400, 0, 4000);
  TH1D C_Up_h_diffangle_track_shower("C_Up_h_diffangle_track_shower", ";diffangle_track_shower", 400, 0, 4000);
  TH1D C_Up_h_myratio50_cascmuon_over_mu("C_Up_h_myratio50_cascmuon_over_mu", ";myratio50_cascmuon_over_mu",400, 0, 1);
  TH1D C_Up_h_myratio30_cascmuon_over_mu("C_Up_h_myratio30_cascmuon_over_mu", ";myratio30_cascmuon_over_mu",400, 0, 1);
  TH1D C_Up_h_myratio50_casc_over_mu("C_Up_h_myratio50_casc_over_mu", ";myratio50_casc_over_mu",400, 0, 1);
  TH1D C_Up_h_myratio30_casc_over_mu("C_Up_h_myratio30_casc_over_mu", ";myratio30_casc_over_mu",400, 0, 1);
  TH1D C_Up_h_ratio_closehits_cascmuon_over_mu("C_Up_h_ratio_closehits_cascmuon_over_mu", ";ratio_closehits_cascmuon_over_mu",400, 0, 1);
  TH1D C_Up_h_ratio_closehits_casc_over_mu("C_Up_h_ratio_closehits_casc_over_mu", ";ratio_closehits_casc_over_mu",400, 0, 1);
  TH1D C_Up_h_redToT_cascmuon_over_mu("C_Up_h_redToT_cascmuon_over_mu", ";redToT_cascmuon_over_mu",400, 0, 1);
  TH1D C_Up_h_redToT_casc_over_mu("C_Up_h_redToT_casc_over_mu", ";redToT_casc_over_mu",400, 0, 1);
  TH1D C_Up_h_bjorken_y("C_Up_h_bjorken_y", ";bjorken_y",400, 0, 1);
  TH1D C_Up_h_diffangle("C_Up_h_diffangle", ";diffangle",180, 0, 180);
  TH1D C_Up_h_LogE_mu ("C_Up_h_LogE_mu","; LogE_mu ",100,0,10);
  TH1D C_Up_h_LogEbundle ("C_Up_h_LogEbundle","; LogEbundle ",100,0,10);
  //TH1DC_Up_//h_logE_mu_max ("//C_Up_//h_logE_mu_max","; LogE_mu_max ",100,0,10);
  TH1D C_Up_h_cos_zen_mu ("C_Up_h_cos_zen_mu","; Cos(zenitC_Up_h_mu)",100,-1,1);
  
  TH2D C_Up_Lik_vs_NNhits ("C_Up_Lik_vs_NNhits","; ",1000,0,1000,700,-100,600);
  TH2D C_Up_Lik_vs_Q1value ("C_Up_Lik_vs_Q1value","; ",500,-5,5,700,-100,600);
  TH2D C_Up_NNhits_vs_Q1value ("C_Up_NNhits_vs_Q1value","; ",500,-5,5,1000,0,1000);
  TH2D C_Up_Nlines_vs_zenith ("C_Up_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D C_Up_Chere_Nlines_vs_zenith ("C_Up_Chere_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D C_Up_TrLengthIT_2_vs_Slen ("C_Up_TrLengthIT_2_vs_Slen",";Slen;TrLengthIT_2",500,0,1000,500,0,1000);
  TH2D C_Up_Diffangle_vs_run ("C_Up_Diffangle_vs_run",";diffangle;run",4667,9333,14000,150,0,150);
  TH2D C_Up_LogEreco_vs_CosZen ("C_Up_LogEreco_vs_CosZen",";CosZen;LogEreco",100,-1,1,100,0,9);
  TH2D C_Up_LogEreco_vs_LogEreco2 ("C_Up_LogEreco_vs_LogEreco2",";LogEreco2;LogEreco",100,-1,10,100,-1,10);
  TH2D C_Up_LogEreco_vs_LogE_mu ("C_Up_LogEreco_vs_LogE_mu",";LogE_mu;LogEreco",100,-1,10,100,-1,10);
  TH2D C_Up_LogEreco2_vs_LogE_mu ("C_Up_LogEreco2_vs_LogE_mu",";LogE_mu;LogEreco2",100,-1,10,100,-1,10);
  TH2D C_Up_Rate3Dmuon ("C_Up_Rate3Dmuon","; run_nr",4667,9333,14000,10000,0,100);
  TH2D C_Up_Rate3Dshower ("C_Up_Rate3Dshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D C_Up_RateMXshower ("C_Up_RateMXshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D C_Up_LogEreco_vs_LogE_neu ("C_Up_LogEreco_vs_LogE_neu",";LogE_neu;LogEreco",100,-1,10,100,-1,10);
  TH2D C_Up_LogEreco2_vs_LogE_neu ("C_Up_LogEreco2_vs_LogE_neu",";LogE_neu;LogEreco2",100,-1,10,100,-1,10);
  TH2D C_Up_LogEresolution_mu_vs_CosZen ("C_Up_LogEresolution_mu_vs_CosZen",";CosZen;(LogEreco-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D C_Up_LogEresolution_mu_cor_vs_CosZen ("C_Up_LogEresolution_mu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D C_Up_LogEresolution_neu_vs_CosZen ("C_Up_LogEresolution_neu_vs_CosZen",";CosZen;(LogEreco-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  TH2D C_Up_LogEresolution_neu_cor_vs_CosZen ("C_Up_LogEresolution_neu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);

  //First cuts
  TH1D First_Up_CosZen ("First_Up_CosZen","; Cos(reco_zenith)",100,-1,1);
  TH1D First_Up_Lik ("First_Up_Lik","; Likelihood", 700,-100,600);
  TH1D First_Up_LogBeta0 ("First_Up_LogBeta0","; LogBeta0",400,-4,0);
  TH1D First_Up_Beta0_deg ("First_Up_Beta0_deg","; LogBeta0deg",600,-4,2);
  TH1D First_Up_DOMs ("First_Up_DOMs","First_Up_DOMs; # of DOMs",200,0,200);
  TH1D First_Up_Cherenkov_DOMs ("First_Up_Cherenkov_DOMs","First_Up_Cherenkov_DOMs; # of Cherenkov DOMs",200,0,200);
  TH1D First_Up_Q1value ("First_Up_Q1value","First_Up_Q1value; -Likelihood/Nhits",500,-5,5);
  TH1D First_Up_NNhits ("First_Up_NNhits","First_Up_NNhits; # of Nhits",1000,0,1000); 
  TH1D First_Up_LogEreco ("First_Up_LogEreco","; LogEreco ",100,0,8);
  TH1D First_Up_TrLen  ("First_Up_TrLen","; Track Length",1000,0,1000);
  TH1D First_Up_Npe  ("First_Up_Npe","; # of p.e",10000,0,10000);
  TH1D First_Up_Ratiopmts  ("First_Up_Ratiopmts","; Chere.PMTs/Reco.PMTs",100,0,1); 
  TH1D First_Up_Ratiodoms  ("First_Up_Ratiodoms","; Chere.DOMs/Reco.DOMs",100,0,1);
  TH2D First_Up_Allreco_rate ("First_Up_Allreco_rate","; run_nr",4667,9333,14000,10000,0,100);
  TH1D First_Up_Rvalue ("First_Up_Rvalue","Up_R; R(m)",1000,0,1000);
  TH1D First_Up_All_triggered_hits ("First_Up_All_triggered_hits","First_Up_All_triggered_hits; # of hits",1000,0,1000);
  TH1D First_Up_Phi ("First_Up_Phi","First_Up_Phi; Azimuth [deg]",360,-180,180);
  TH1D First_Up_Zenith_shower ("First_Up_Zenith_shower","First_Up_zenith_shower; [deg]",360,0,360);
  TH1D First_Up_DLik ("First_Up_DLik","; First_Up_DLikelihood", 700,-100,600);
  TH1D First_Up_Ratio6 ("First_Up_Ratio6","; First_Up_Ratio6",1200,0,1.2);
  TH1D First_Up_MaxZenSol ("First_Up_MaxZenSol","First_Up_MaxZenSol; [deg]",360,0,360);
  TH1D First_Up_RedToT_muon ("First_Up_RedToT_muon","; Total_tot/nhits",300,0,300);
  TH1D First_Up_Zenith ("First_Up_Zenith","First_Up_Zenith; [deg]",360,0,360);
  TH1D First_Up_Delta_ZeniTh_track_shower ("First_Up_Delta_ZeniTh_track_shower","First_Up_Delta_ZeniTh_track_shower; [deg]",360,-180,180);
  TH1D First_Up_MinZenSol ("First_Up_MinZenSol","First_Up_MinZenSol; [deg]",360,0,360);
  TH1D First_Up_ItoLen ("First_Up_ItoLen","; First_Up_ItoLen",1200,0,1.2);
  TH1D First_Up_Ratio330 ("First_Up_Ratio330","; First_Up_Ratio330",1200,0,1.2);
  TH1D First_Up_Num_good_sol ("First_Up_Num_good_sol",";# of good solutions",100,0,100);
  TH1D First_Up_PosZ  ("First_Up_PosZ","; First_Up_PosZ",1000,0,1000);
  TH1D First_Up_LogEnu ("First_Up_LogEnu","; First_Up_LogEnu ",100,0,8);
  TH1D First_Up_LogAEnergy ("First_Up_LogAEnergy","; First_Up_LogAEnergy ",100,0,8);
  TH1D First_Up_NLines ("First_Up_NLines","; # of lines",12,0,12);
  TH1D First_Up_LogEreco_cor ("First_Up_LogEreco_cor","; LogEreco_cor ",100,0,8);
  TH1D First_Up_Chere_NLines ("First_Up_Chere_NLines","; # of lines",12,0,12);
  //TH1D First_Up_Reco_Dclosest ("First_Up_Reco_Dclosest","; d_closest ", 1000, 0, 1000);
  //TH1D First_Up_Reco_tres ("First_Up_Reco_tres","; reco_tres ", 1000, -500, 500);
  //TH1D First_Up_Reco_Cos_Angle ("First_Up_Reco_Cos_Angle","; reco_cos_angle ", 500, -360, 360);
  TH1D First_Up_Ratio_upper_hemisphere ("First_Up_Ratio_upper_hemisphere","; upper_hits/all_hits",100,0,1);
  TH1D First_Up_Ratio_lower_hemisphere ("First_Up_Ratio_lower_hemisphere","; lowerr_hits/all_hits",100,0,1);
  TH1D First_Up_Ratio_border_upper_hemisphere ("First_Up_Ratio_border_upper_hemisphere","; border_upper_hits/all_border_hits",100,0,1);
  TH1D First_Up_Ratio_border_lower_hemisphere ("First_Up_Ratio_border_lower_hemisphere","; border_lowerr_hits/all_border_hits",100,0,1);
  TH1D First_Up_Ratio1 ("First_Up_Ratio1","; Ratio1",1200,0,1.2);
  TH1D First_Up_Ratio110 ("First_Up_Ratio110","; Ratio110",1200,0,1.2);
  TH1D First_Up_Ratio130 ("First_Up_Ratio130","; Ratio130",1200,0,1.2);
  TH1D First_Up_NtrackIT_ratio ("First_Up_NtrackIT_ratio","; NtrackIT/Nhits",1200,0,1.2);
  TH1D First_Up_TrLenIT_3  ("First_Up_TrLenIT_3","; Track Length",1000,0,1000);

  TH1D  First_Up_h_ToT_border_mu("First_Up_h_ToT_border_mu",";ToT_border_mu ",400,0,4000);
  TH1D  First_Up_h_ToT_border_casc("First_Up_h_ToT_border_casc",";ToT_border_casc ",400,0,4000);
  TH1D  First_Up_h_ToT_border_cascmu("First_Up_h_ToT_border_cascmu",";ToT_border_cascmu ",400,0,4000);
  TH1D  First_Up_h_ToT_trig("First_Up_h_ToT_trig",";ToT_trig ",400,0,4000);
  TH1D  First_Up_h_max_ToT_trig("First_Up_h_max_ToT_trig",";max_ToT_trig ",400,0,4000);
  TH1D  First_Up_h_ToT_IT("First_Up_h_ToT_IT",";ToT_IT ",400,0,4000);
  TH1D  First_Up_h_ToT_allIT("First_Up_h_ToT_allIT",";ToT_allIT ",400,0,4000);
  TH1D  First_Up_h_sum_ToT_casc("First_Up_h_sum_ToT_casc",";sum_ToT_casc" ,400,0,4000);
  TH1D  First_Up_h_sum_ToT_mu("First_Up_h_sum_ToT_mu",";sum_ToT_mu" ,400,0,4000);
  TH1D  First_Up_h_sum_ToT_casc_mu("First_Up_h_sum_ToT_casc_mu",";sum_ToT_casc_mu" ,400,0,4000);
  //----------------------------------------------------------
  TH1D First_Up_h_Ratio430 ("First_Up_h_Ratio430","; Ratio430",1200,0,1.2);
  TH1D First_Up_h_Ratio5 ("First_Up_h_Ratio5","; Ratio5",1200,0,1.2);
  TH1D First_Up_h_DiffDistance_mu  ("First_Up_h_DiffDistance_mu","; dmax-dmin [m]",1400,0,1400);
  TH1D First_Up_h_RatioCloseHits_mu  ("First_Up_h_RatioCloseHits_mu","; nhits_400m/nhits",400,0,1);
  TH1D First_Up_h_NtrackEarly("First_Up_h_NtrackEarly", ";NtrackEarly",400,0,4000);
  TH1D First_Up_h_D("First_Up_h_D", ";D", 400, 0, 4000);
  TH1D First_Up_h_NtrackIT30("First_Up_h_NtrackIT30", ";NtrackIT30", 400, 0, 4000);
  TH1D First_Up_h_NtrackLate("First_Up_h_NtrackLate", ";NtrackLate", 400, 0, 4000);
  TH1D First_Up_h_myratio50_muon("First_Up_h_myratio50_muon", ";myratio50_muon", 400, 0, 1);
  TH1D First_Up_h_ratio_cherenkov_lines("First_Up_h_ratio_cherenkov_lines", ";ratio_cherenkov_lines", 400, 0, 1);
  TH1D First_Up_h_SnpeT("First_Up_h_SnpeT", ";SnpeT", 400, 0, 4000);
  TH1D First_Up_h_mean_tres_it("First_Up_h_mean_tres_it", ";mean_tres_it", 400, 0, 4000);
  TH1D First_Up_h_max_lik_up("First_Up_h_max_lik_up", ";max_lik_up", 400, 0, 4000);
  TH1D First_Up_h_max_lik_down("First_Up_h_max_lik_down", ";max_lik_down", 400, 0, 4000);
  TH1D First_Up_h_diff_theta("First_Up_h_diff_theta", ";diff_theta", 400, 0, 4000);
  TH1D First_Up_h_diff_dist_casc_mu("First_Up_h_diff_dist_casc_mu", ";diff_dist_casc_mu", 400, 0, 4000);
  TH1D First_Up_h_diff_dist_casc("First_Up_h_diff_dist_casc", ";diff_dist_casc", 400, 0, 4000);
  TH1D First_Up_h_ratio_closehits_cascmuon("First_Up_h_ratio_closehits_cascmuon", ";ratio_closehits_cascmuon", 400, 0, 1);
  TH1D First_Up_h_ratio_closehits_casc("First_Up_h_ratio_closehits_casc", ";ratio_closehits_casc", 400, 0,1);
  TH1D First_Up_h_redToT_cascmuon("First_Up_h_redToT_cascmuon", ";redToT_cascmuon", 400, 0, 4000);
  TH1D First_Up_h_redToT_casc("First_Up_h_redToT_casc", ";redToT_casc", 400, 0, 4000);
  TH1D First_Up_h_myratio50_cascmuon("First_Up_h_myratio50_cascmuon", ";myratio50_cascmuon", 400, 0, 1);
  TH1D First_Up_h_myratio50_casc("First_Up_h_myratio50_casc", ";myratio50_casc", 400, 0, 1);
  TH1D First_Up_h_myratio30_cascmuon("First_Up_h_myratio30_cascmuon", ";myratio30_cascmuon", 400, 0, 1);
  TH1D First_Up_h_myratio30_muon("First_Up_h_myratio30_muon", ";myratio30_muon", 400, 0, 1);
  TH1D First_Up_h_myratio30_casc("First_Up_h_myratio30_casc", ";myratio30_casc", 400, 0, 1);
  TH1D First_Up_h_min_diff_sollik("First_Up_h_min_diff_sollik",";min_diff_sollik", 400, 0, 4000);
  TH1D First_Up_h_beta0_shower_deg("First_Up_h_beta0_shower_deg", ";beta0_shower_deg", 400, 0, 4000);
  TH1D First_Up_h_lik_shower("First_Up_h_lik_shower", ";lik_shower", 400, 0, 4000);
  TH1D First_Up_h_best_trk_pos_shower_z("First_Up_h_best_trk_pos_shower_z", ";best_trk_pos_shower_z", 400, 0, 4000);
  TH1D First_Up_h_normdlik("First_Up_h_normdlik", ";normdlik", 400, 0, 4000);
  TH1D First_Up_h_upsol("First_Up_h_upsol", ";upsol", 100, 0, 60);
  TH1D First_Up_h_Slen("First_Up_h_Slen", ";Slen", 1000, 0, 1000);
  TH1D First_Up_h_max_diff_sollik("First_Up_h_max_diff_sollik", ";max_diff_sollik", 400, 0, 4000);
  TH1D First_Up_h_diffangle_track_shower("First_Up_h_diffangle_track_shower", ";diffangle_track_shower", 400, 0, 4000);
  TH1D First_Up_h_myratio50_cascmuon_over_mu("First_Up_h_myratio50_cascmuon_over_mu", ";myratio50_cascmuon_over_mu",400, 0, 1);
  TH1D First_Up_h_myratio30_cascmuon_over_mu("First_Up_h_myratio30_cascmuon_over_mu", ";myratio30_cascmuon_over_mu",400, 0, 1);
  TH1D First_Up_h_myratio50_casc_over_mu("First_Up_h_myratio50_casc_over_mu", ";myratio50_casc_over_mu",400, 0, 1);
  TH1D First_Up_h_myratio30_casc_over_mu("First_Up_h_myratio30_casc_over_mu", ";myratio30_casc_over_mu",400, 0, 1);
  TH1D First_Up_h_ratio_closehits_cascmuon_over_mu("First_Up_h_ratio_closehits_cascmuon_over_mu", ";ratio_closehits_cascmuon_over_mu",400, 0, 1);
  TH1D First_Up_h_ratio_closehits_casc_over_mu("First_Up_h_ratio_closehits_casc_over_mu", ";ratio_closehits_casc_over_mu",400, 0, 1);
  TH1D First_Up_h_redToT_cascmuon_over_mu("First_Up_h_redToT_cascmuon_over_mu", ";redToT_cascmuon_over_mu",400, 0, 1);
  TH1D First_Up_h_redToT_casc_over_mu("First_Up_h_redToT_casc_over_mu", ";redToT_casc_over_mu",400, 0, 1);
  TH1D First_Up_h_bjorken_y("First_Up_h_bjorken_y", ";bjorken_y",400, 0, 1);
  TH1D First_Up_h_diffangle("First_Up_h_diffangle", ";diffangle",180, 0, 180);
  TH1D First_Up_h_LogE_mu ("First_Up_h_LogE_mu","; LogE_mu ",100,0,10);
  TH1D First_Up_h_LogEbundle ("First_Up_h_LogEbundle","; LogEbundle ",100,0,10);
  //TH1DFirst_Up_//h_logE_mu_max ("//First_Up_//h_logE_mu_max","; LogE_mu_max ",100,0,10);
  TH1D First_Up_h_cos_zen_mu ("First_Up_h_cos_zen_mu","; Cos(zenitFirst_Up_h_mu)",100,-1,1);
 
  TH2D First_Up_Lik_vs_NNhits ("First_Up_Lik_vs_NNhits","; ",1000,0,1000,700,-100,600);
  TH2D First_Up_Lik_vs_Q1value ("First_Up_Lik_vs_Q1value","; ",500,-5,5,700,-100,600);
  TH2D First_Up_NNhits_vs_Q1value ("First_Up_NNhits_vs_Q1value","; ",500,-5,5,1000,0,1000);
  TH2D First_Up_Nlines_vs_zenith ("First_Up_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D First_Up_Chere_Nlines_vs_zenith ("First_Up_Chere_Nlines_vs_zenith","; ",360,0,360,12,0,12);
  TH2D First_Up_TrLengthIT_2_vs_Slen ("First_Up_TrLengthIT_2_vs_Slen",";Slen;TrLengthIT_2",500,0,1000,500,0,1000);
  TH2D First_Up_Diffangle_vs_run ("First_Up_Diffangle_vs_run",";diffangle;run",4667,9333,14000,150,0,150);
  TH2D First_Up_LogEreco_vs_CosZen ("First_Up_LogEreco_vs_CosZen",";CosZen;LogEreco",100,-1,1,100,0,9);
  TH2D First_Up_LogEreco_vs_LogEreco2 ("First_Up_LogEreco_vs_LogEreco2",";LogEreco2;LogEreco",100,-1,10,100,-1,10);
  TH2D First_Up_LogEreco_vs_LogE_mu ("First_Up_LogEreco_vs_LogE_mu",";LogE_mu;LogEreco",100,-1,10,100,-1,10);
  TH2D First_Up_LogEreco2_vs_LogE_mu ("First_Up_LogEreco2_vs_LogE_mu",";LogE_mu;LogEreco2",100,-1,10,100,-1,10);
  TH2D First_Up_Rate3Dmuon ("First_Up_Rate3Dmuon","; run_nr",4667,9333,14000,10000,0,100);
  TH2D First_Up_Rate3Dshower ("First_Up_Rate3Dshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D First_Up_RateMXshower ("First_Up_RateMXshower","; run_nr",4667,9333,14000,10000,0,100);
  TH2D First_Up_LogEreco_vs_LogE_neu ("First_Up_LogEreco_vs_LogE_neu",";LogE_neu;LogEreco",100,-1,10,100,-1,10);
  TH2D First_Up_LogEreco2_vs_LogE_neu ("First_Up_LogEreco2_vs_LogE_neu",";LogE_neu;LogEreco2",100,-1,10,100,-1,10);
  TH2D First_Up_LogEresolution_mu_vs_CosZen ("First_Up_LogEresolution_mu_vs_CosZen",";CosZen;(LogEreco-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D First_Up_LogEresolution_mu_cor_vs_CosZen ("First_Up_LogEresolution_mu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D First_Up_LogEresolution_neu_vs_CosZen ("First_Up_LogEresolution_neu_vs_CosZen",";CosZen;(LogEreco-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);
  TH2D First_Up_LogEresolution_neu_cor_vs_CosZen ("First_Up_LogEresolution_neu_cor_vs_CosZen",";CosZen;(LogEreco_cor-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);

  //Common level histos
  TH1D Common_CosZen ("Common_CosZen","; Cos(reco_zenith)",100,-1,1);     
  TH1D Common_Lik ("Common_Lik","; Likelihood", 700,-100,600);
  TH1D Common_LogBeta0 ("Common_LogBeta0","; LogBeta0",400,-4,0);
  TH1D Common_LogEreco ("Common_LogEreco","; LogEreco ",100,0,8);
  TH1D Common_TrLen  ("Common_TrLen","; Track Length",1000,0,1000);
  TH2D Common_LogEresolution_mu_vs_CosZen ("Common_LogEresolution_mu_vs_CosZen",";CosZen;(LogEreco-LogE_mu) / LogE_mu",100,-1,1,100,-3,3);
  TH2D Common_LogEresolution_neu_vs_CosZen ("Common_LogEresolution_neu_vs_CosZen",";CosZen;(LogEreco-LogE_neu) / LogE_neu",100,-1,1,100,-3,3);

  
  //==========================================================================================================
  //==========================================================================================================

  //=================================================================
  //======== Start the analysis =====================================
  //=================================================================
  int current_run_number=0;
  int all_evts,pass_evts,pass_up_evts,pass_C_evts,pass_C_up_evts,pass_First_evts,pass_First_up_evts;
  all_evts=pass_evts=pass_up_evts=pass_C_evts=pass_C_up_evts=pass_First_evts=pass_First_up_evts=0;
  float reco_events,up_reco_events,c_reco_events,c_up_reco_events,first_reco_events,first_up_reco_events;
  float Muon_evts,up_Muon_evts,c_Muon_evts,c_up_Muon_evts,first_Muon_evts,first_up_Muon_evts;
  float Shower_evts,up_Shower_evts,c_Shower_evts,c_up_Shower_evts,first_Shower_evts,first_up_Shower_evts;
  float MXShower_evts,up_MXShower_evts,c_MXShower_evts,c_up_MXShower_evts,first_MXShower_evts,first_up_MXShower_evts;
  
  //loop over all entries (each entry is an event NOT a run) (ProcessedEvents tree)
  for( int i=0 ; i<dataTree->GetEntries() ; i++ ){


    all_evts++;

    dataTree->GetEntry(i); // get i-th event parameters

    //============================= new variables ====================================//
    float logbeta0_deg;
    logbeta0_deg=TMath::Log10(jbeta0_deg);
     
    float Nhits;
    Nhits = GNhit;

    float chere_nlines;
    chere_nlines = static_cast<float>(num_cherenkov_lines);

    float_num_triggered_doms=static_cast<float>(num_triggered_doms);

    double logEresolution_mu, logEresolution_cor_mu, logEresolution_neu, logEresolution_cor_neu;
    logEresolution_mu =( logEreco - logE_mu ) / logE_mu;
    logEresolution_cor_mu = ( logEreco2 - logE_mu ) / logE_mu;
    logEresolution_neu =( logEreco - logE_nu ) / logE_nu;
    logEresolution_cor_neu = ( logEreco2 - logE_nu ) / logE_nu;
    //====================================================================================//

    //Define cuts. In order to work properly, for each loop set them false again
    bool anti_noise_cuts=false;
    if(jlik > 50 && Nhits > 20){anti_noise_cuts=true;}

    bool first_precuts=false;
    if(jlik > 40 && Slen > 100 && zenith > 80 && logbeta0 <-1.5){
    if(logEreco > 2.7 && jlik >50.0 && TrLengthIT_2 >100.0 && logbeta0 <-1.5 && num_triggered_doms>=5){
      first_precuts=true;}}

    bool common_level=false;
    if(jlik > 40 && Slen > 100 && zenith > 80 && logbeta0 <-1.5){common_level=true;}

    bool up_cut=false;
    if(zenith > 90){up_cut=true;}

    //extra condition for diffangle WR <5 br >10
    if(1==1){//diffangle < 5){
     
      //==== how many reco events for each run  ========//
      //try to distinguish each event to its respective run
      if(pseudo_runid != current_run_number){
	current_run_number=pseudo_runid;
	reco_events=up_reco_events=c_reco_events=c_up_reco_events=first_reco_events=first_up_reco_events=0;
        Muon_evts=up_Muon_evts=c_Muon_evts=c_up_Muon_evts=first_Muon_evts=first_up_Muon_evts=0;
	Shower_evts=up_Shower_evts=c_Shower_evts=c_up_Shower_evts=first_Shower_evts=first_up_Shower_evts=0;
	MXShower_evts=up_MXShower_evts=c_MXShower_evts=c_up_MXShower_evts=first_MXShower_evts=first_up_MXShower_evts=0;
      }

      
      reco_events++;
      if(flag_muon_3D > 0){Muon_evts++;}
      if(flag_shower_3D > 0){Shower_evts++;}
      if(flag_shower_MX > 0){MXShower_evts++;}
      reco_events_to_run[pseudo_runid]=reco_events;
      evts_3d_muon_to_run[pseudo_runid]=Muon_evts;
      evts_3d_shower_to_run[pseudo_runid]=Shower_evts;
      evts_mx_shower_to_run[pseudo_runid]=MXShower_evts;
      pass_evts++;//

      //std::cout<<"3dmuon: "<<flag_muon_3D<<" show: "<<flag_shower_3D<<" mx: "<<flag_shower_MX<<std::endl;
      
      if(up_cut){up_reco_events++;
	if(flag_muon_3D){up_Muon_evts++;}
	if(flag_shower_3D){up_Shower_evts++;}
	if(flag_shower_MX){up_MXShower_evts++;}
	up_reco_events_to_run[pseudo_runid]=up_reco_events;
	up_evts_3d_muon_to_run[pseudo_runid]=up_Muon_evts;
	up_evts_3d_shower_to_run[pseudo_runid]=up_Shower_evts;
	up_evts_mx_shower_to_run[pseudo_runid]=up_MXShower_evts;
	pass_up_evts++;//
      }
   
      if(anti_noise_cuts){
	c_reco_events++;
	if(flag_muon_3D){c_Muon_evts++;}
	if(flag_shower_3D){c_Shower_evts++;}
	if(flag_shower_MX){c_MXShower_evts++;}
	c_reco_events_to_run[pseudo_runid]=c_reco_events;
	c_evts_3d_muon_to_run[pseudo_runid]=c_Muon_evts;
	c_evts_3d_shower_to_run[pseudo_runid]=c_Shower_evts;
	c_evts_mx_shower_to_run[pseudo_runid]=c_MXShower_evts;
	pass_C_evts++;//
	if(up_cut){c_up_reco_events++;
	  if(flag_muon_3D){c_up_Muon_evts++;}
	  if(flag_shower_3D){c_up_Shower_evts++;}
	  if(flag_shower_MX){c_up_MXShower_evts++;}
	  c_up_reco_events_to_run[pseudo_runid]=c_up_reco_events;
	  c_up_evts_3d_muon_to_run[pseudo_runid]=c_up_Muon_evts;
	  c_up_evts_3d_shower_to_run[pseudo_runid]=c_up_Shower_evts;
	  c_up_evts_mx_shower_to_run[pseudo_runid]=c_up_MXShower_evts;
	  pass_C_up_evts++;//
	}
      }
      
      if(first_precuts){
	first_reco_events++;
	if(flag_muon_3D){first_Muon_evts++;}
	if(flag_shower_3D){first_Shower_evts++;}
	if(flag_shower_MX){first_MXShower_evts++;}
	first_reco_events_to_run[pseudo_runid]=first_reco_events;
	first_evts_3d_muon_to_run[pseudo_runid]=first_Muon_evts;
	first_evts_3d_shower_to_run[pseudo_runid]=first_Shower_evts;
	first_evts_mx_shower_to_run[pseudo_runid]=first_MXShower_evts;
	pass_First_evts++;//
	if(up_cut){first_up_reco_events++;
	  if(flag_muon_3D){first_up_Muon_evts++;}
	  if(flag_shower_3D){first_up_Shower_evts++;}
	  if(flag_shower_MX){first_up_MXShower_evts++;}
	  first_up_reco_events_to_run[pseudo_runid]=first_up_reco_events;
	  first_up_evts_3d_muon_to_run[pseudo_runid]=first_up_Muon_evts;
	  first_up_evts_3d_shower_to_run[pseudo_runid]=first_up_Shower_evts;
	  first_up_evts_mx_shower_to_run[pseudo_runid]=first_up_MXShower_evts;
	  pass_First_up_evts++;//
	}
      }
      //==========================================================//
        
      double the_w2 = w3ice / icecube_flux_diff2019(E_nu);
      double new_w3ice = the_w2 * icecube_flux_diff2021(E_nu);
    				
      //double w3 = 1.0;
      double w3 = new_w3ice;
      //double w3 = w3atm;
      //--- Every value has anti-bug cuts from the pre-process
      //std::cout<<pseudo_runid<<std::endl;
     
      CosZen.Fill(cos_zen,w3);
      Lik.Fill(jlik,w3);
      LogBeta0.Fill(logbeta0,w3);
      Beta0_deg.Fill(logbeta0_deg,w3);
      DOMs.Fill(num_triggered_doms,w3);
      Cherenkov_DOMs.Fill(num_cherenkov_doms,w3);
      Q1value.Fill(Q1,w3);
      NNhits.Fill(Nhits,w3);
      LogEreco.Fill(logEreco,w3);    
      TrLen.Fill(TrLengthIT_2,w3);
      Npe.Fill(SnpeT,w3);
      Ratiopmts.Fill(ratio_cherenkov_pmts,w3);
      Ratiodoms.Fill(ratio_cherenkov_doms,w3);
      Rvalue.Fill(R,w3);
      All_triggered_hits.Fill(num_triggered_hits,w3);
      Phi.Fill(fi,w3);
      Zenith_shower.Fill(zenith_shower,w3);
      DLik.Fill(dlik,w3);
      Ratio6.Fill(ratio6,w3);
      MaxZenSol.Fill(max_zen_sol,w3);
      RedToT_muon.Fill(redToT_muon,w3);
      Zenith.Fill(zenith,w3);
      Delta_ZeniTh_track_shower.Fill(delta_zenith_track_shower,w3);
      MinZenSol.Fill(min_zen_sol,w3);
      ItoLen.Fill(itoverlen,w3);
      Ratio330.Fill(ratio330,w3);
      Num_good_sol.Fill(num_of_good_sol,w3);
      PosZ.Fill(best_trk_pos_z,w3);
      LogEnu.Fill(logE_nu,w3);
      LogAEnergy.Fill(logEreco_shower,w3);
      NLines.Fill(num_triggered_lines,w3);
      LogEreco_cor.Fill(logEreco2,w3);
      Chere_NLines.Fill(num_cherenkov_lines,w3);
      Ratio_upper_hemisphere.Fill(Nhits_upper/(Nhits_upper + Nhits_lower), w3);
      Ratio_lower_hemisphere.Fill(Nhits_lower/(Nhits_upper + Nhits_lower), w3);
      Ratio_border_upper_hemisphere.Fill(Nhits_border_upper/Nborder_hits, w3);
      Ratio_border_lower_hemisphere.Fill(Nhits_border_lower/Nborder_hits, w3);
      Ratio1.Fill(ratio1,w3);
      Ratio110.Fill(ratio110,w3);
      Ratio130.Fill(ratio130,w3);
      NtrackIT_ratio.Fill(NtrackIT/Nhits, w3); //number of in time hits over all
      TrLenIT_3.Fill(TrLengthIT_3, w3);

      h_ToT_border_mu.Fill(ToT_border_mu,w3);
      h_ToT_border_casc.Fill(ToT_border_casc,w3);
      h_ToT_border_cascmu.Fill(ToT_border_cascmu,w3);
      h_ToT_trig.Fill(ToT_trig,w3);
      h_max_ToT_trig.Fill(max_ToT_trig,w3);
      h_ToT_IT.Fill(ToT_IT,w3);
      h_ToT_allIT.Fill(ToT_allIT,w3);
      h_sum_ToT_casc.Fill(sum_ToT_casc,w3);
      h_sum_ToT_mu.Fill(sum_ToT_mu,w3);
      h_sum_ToT_casc_mu.Fill(sum_ToT_casc_mu,w3);

      //-----------------------------------------------
      h_Ratio430.Fill(ratio430,w3);
      h_Ratio5.Fill(ratio5,w3);
      h_DiffDistance_mu.Fill(diff_dist_mu,w3);
      h_RatioCloseHits_mu.Fill(ratio_closehits_muon,w3);
      h_NtrackEarly.Fill(NtrackEarly,w3);
      h_D.Fill(D,w3);
      h_NtrackIT30.Fill(NtrackIT30,w3);
      h_NtrackLate.Fill(NtrackLate,w3);
      h_myratio50_muon.Fill(myratio50_muon,w3);
      h_ratio_cherenkov_lines.Fill(ratio_cherenkov_lines,w3);
      h_SnpeT.Fill(SnpeT,w3);
      h_mean_tres_it.Fill(mean_tres_it,w3);
      h_max_lik_up.Fill(max_lik_up,w3);
      h_max_lik_down.Fill(max_lik_down,w3);
      h_diff_theta.Fill(diff_theta,w3);
      h_diff_dist_casc_mu.Fill(diff_dist_casc_mu,w3);
      h_diff_dist_casc.Fill(diff_dist_casc,w3);
      h_ratio_closehits_cascmuon.Fill(ratio_closehits_cascmuon,w3);
      h_ratio_closehits_casc.Fill(ratio_closehits_casc,w3);
      h_redToT_cascmuon.Fill(redToT_cascmuon,w3);
      h_redToT_casc.Fill(redToT_casc,w3);
      h_myratio50_cascmuon.Fill(myratio50_cascmuon,w3);
      h_myratio50_casc.Fill(myratio50_casc,w3);
      h_myratio30_cascmuon.Fill(myratio30_cascmuon,w3);
      h_myratio30_muon.Fill(myratio30_muon,w3);
      h_myratio30_casc.Fill(myratio30_casc,w3);
      h_min_diff_sollik.Fill(min_diff_sollik,w3);
      h_beta0_shower_deg.Fill(beta0_shower_deg,w3);
      h_lik_shower.Fill(lik_shower,w3);
      h_best_trk_pos_shower_z.Fill(best_trk_pos_shower_z,w3);
      h_normdlik.Fill(normdlik,w3);
      h_upsol.Fill(upsol,w3);
      h_Slen.Fill(Slen,w3);
      h_max_diff_sollik.Fill(max_diff_sollik,w3);
      h_diffangle_track_shower.Fill(diffangle_track_shower,w3);
      h_myratio50_cascmuon_over_mu.Fill(myratio50_cascmuon_over_mu,w3);
      h_myratio30_cascmuon_over_mu.Fill(myratio30_cascmuon_over_mu,w3);
      h_myratio50_casc_over_mu.Fill(myratio50_casc_over_mu,w3);
      h_myratio30_casc_over_mu.Fill(myratio30_casc_over_mu,w3);
      h_ratio_closehits_cascmuon_over_mu.Fill(ratio_closehits_cascmuon_over_mu,w3);
      h_ratio_closehits_casc_over_mu.Fill(ratio_closehits_casc_over_mu,w3);
      h_redToT_cascmuon_over_mu.Fill(redToT_cascmuon_over_mu,w3);
      h_redToT_casc_over_mu.Fill(redToT_casc_over_mu,w3);
      h_bjorken_y.Fill(bjorken_y,w3);
      h_diffangle.Fill(diffangle,w3);
      h_LogE_mu.Fill(logE_mu,w3);
      h_LogEbundle.Fill(logEbundle,w3);
      //h_logE_mu_max.Fill(TMath::Log10(E_mu_max),w3);
      h_cos_zen_mu.Fill(cos_zen_mu,w3);
      //----------------------------------------------


      Lik_vs_NNhits.Fill(Nhits,jlik,w3);
      Lik_vs_Q1value.Fill(Q1,jlik,w3);
      NNhits_vs_Q1value.Fill(Q1,Nhits,w3);
      Nlines_vs_zenith.Fill(zenith,num_triggered_lines,w3);
      Chere_Nlines_vs_zenith.Fill(zenith,num_cherenkov_lines,w3);
      TrLengthIT_2_vs_Slen.Fill(Slen,TrLengthIT_2,w3);
      Diffangle_vs_run.Fill(pseudo_runid,diffangle,w3);
      LogEreco_vs_CosZen.Fill(cos_zen,logEreco,w3);
      LogEreco_vs_LogEreco2.Fill(logEreco2,logEreco,w3);
      LogEreco_vs_LogE_mu.Fill(logE_mu,logEreco,w3);
      LogEreco2_vs_LogE_mu.Fill(logE_mu,logEreco2,w3);
      LogEreco_vs_LogE_neu.Fill(logE_nu,logEreco,w3);
      LogEreco2_vs_LogE_neu.Fill(logE_nu,logEreco2,w3);
      LogEresolution_mu_vs_CosZen.Fill(cos_zen,logEresolution_mu,w3);
      LogEresolution_mu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_mu,w3);
      LogEresolution_neu_vs_CosZen.Fill(cos_zen,logEresolution_neu,w3);
      LogEresolution_neu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_neu,w3);
 
      //up
      if(up_cut){	
        Up_CosZen.Fill(cos_zen,w3);
        Up_Lik.Fill(jlik,w3);
        Up_LogBeta0.Fill(logbeta0,w3);
        Up_Beta0_deg.Fill(logbeta0_deg,w3);
        Up_DOMs.Fill(num_triggered_doms,w3);
        Up_Cherenkov_DOMs.Fill(num_cherenkov_doms,w3);
        Up_Q1value.Fill(Q1,w3);
        Up_NNhits.Fill(Nhits,w3);
        Up_LogEreco.Fill(logEreco,w3);
        Up_TrLen.Fill(TrLengthIT_2,w3);
        Up_Npe.Fill(SnpeT,w3);
        Up_Ratiopmts.Fill(ratio_cherenkov_pmts,w3);
        Up_Ratiodoms.Fill(ratio_cherenkov_doms,w3);
        Up_Rvalue.Fill(R,w3);
        Up_All_triggered_hits.Fill(num_triggered_hits,w3);
        Up_Phi.Fill(fi,w3);
	Up_Zenith_shower.Fill(zenith_shower,w3);
        Up_DLik.Fill(dlik,w3);
	Up_Ratio6.Fill(ratio6,w3);
	Up_MaxZenSol.Fill(max_zen_sol,w3);
	Up_RedToT_muon.Fill(redToT_muon,w3);
	Up_Zenith.Fill(zenith,w3);
	Up_Delta_ZeniTh_track_shower.Fill(delta_zenith_track_shower,w3);
	Up_MinZenSol.Fill(min_zen_sol,w3);
	Up_ItoLen.Fill(itoverlen,w3);
	Up_Ratio330.Fill(ratio330,w3);
	Up_Num_good_sol.Fill(num_of_good_sol,w3);
	Up_PosZ.Fill(best_trk_pos_z,w3);
	Up_LogEnu.Fill(logE_nu,w3);
	Up_LogAEnergy.Fill(logEreco_shower,w3);
	Up_NLines.Fill(num_triggered_lines,w3);
	Up_LogEreco_cor.Fill(logEreco2,w3);
	Up_Chere_NLines.Fill(num_cherenkov_lines,w3);
	Up_Ratio_upper_hemisphere.Fill(Nhits_upper/(Nhits_upper + Nhits_lower), w3);
	Up_Ratio_lower_hemisphere.Fill(Nhits_lower/(Nhits_upper + Nhits_lower), w3);
	Up_Ratio_border_upper_hemisphere.Fill(Nhits_border_upper/Nborder_hits, w3);
	Up_Ratio_border_lower_hemisphere.Fill(Nhits_border_lower/Nborder_hits, w3);
	Up_Ratio1.Fill(ratio1,w3);
	Up_Ratio110.Fill(ratio110,w3);
	Up_Ratio130.Fill(ratio130,w3);
	Up_NtrackIT_ratio.Fill(NtrackIT/Nhits, w3);
	Up_TrLenIT_3.Fill(TrLengthIT_3, w3);

	Up_h_ToT_border_mu.Fill(ToT_border_mu,w3);
	Up_h_ToT_border_casc.Fill(ToT_border_casc,w3);
	Up_h_ToT_border_cascmu.Fill(ToT_border_cascmu,w3);
	Up_h_ToT_trig.Fill(ToT_trig,w3);
	Up_h_max_ToT_trig.Fill(max_ToT_trig,w3);
	Up_h_ToT_IT.Fill(ToT_IT,w3);
	Up_h_ToT_allIT.Fill(ToT_allIT,w3);
	Up_h_sum_ToT_casc.Fill(sum_ToT_casc,w3);
	Up_h_sum_ToT_mu.Fill(sum_ToT_mu,w3);
	Up_h_sum_ToT_casc_mu.Fill(sum_ToT_casc_mu,w3);

	//------------------------------------------------
	Up_h_Ratio430.Fill(ratio430,w3);
	Up_h_Ratio5.Fill(ratio5,w3);
	Up_h_DiffDistance_mu.Fill(diff_dist_mu,w3);
	Up_h_RatioCloseHits_mu.Fill(ratio_closehits_muon,w3);
	Up_h_NtrackEarly.Fill(NtrackEarly,w3);
	Up_h_D.Fill(D,w3);
	Up_h_NtrackIT30.Fill(NtrackIT30,w3);
	Up_h_NtrackLate.Fill(NtrackLate,w3);
	Up_h_myratio50_muon.Fill(myratio50_muon,w3);
	Up_h_ratio_cherenkov_lines.Fill(ratio_cherenkov_lines,w3);
	Up_h_SnpeT.Fill(SnpeT,w3);
	Up_h_mean_tres_it.Fill(mean_tres_it,w3);
	Up_h_max_lik_up.Fill(max_lik_up,w3);
	Up_h_max_lik_down.Fill(max_lik_down,w3);
	Up_h_diff_theta.Fill(diff_theta,w3);
	Up_h_diff_dist_casc_mu.Fill(diff_dist_casc_mu,w3);
	Up_h_diff_dist_casc.Fill(diff_dist_casc,w3);
	Up_h_ratio_closehits_cascmuon.Fill(ratio_closehits_cascmuon,w3);
	Up_h_ratio_closehits_casc.Fill(ratio_closehits_casc,w3);
	Up_h_redToT_cascmuon.Fill(redToT_cascmuon,w3);
	Up_h_redToT_casc.Fill(redToT_casc,w3);
	Up_h_myratio50_cascmuon.Fill(myratio50_cascmuon,w3);
	Up_h_myratio50_casc.Fill(myratio50_casc,w3);
	Up_h_myratio30_cascmuon.Fill(myratio30_cascmuon,w3);
	Up_h_myratio30_muon.Fill(myratio30_muon,w3);
	Up_h_myratio30_casc.Fill(myratio30_casc,w3);
	Up_h_min_diff_sollik.Fill(min_diff_sollik,w3);
	Up_h_beta0_shower_deg.Fill(beta0_shower_deg,w3);
	Up_h_lik_shower.Fill(lik_shower,w3);
	Up_h_best_trk_pos_shower_z.Fill(best_trk_pos_shower_z,w3);
	Up_h_normdlik.Fill(normdlik,w3);
	Up_h_upsol.Fill(upsol,w3);
	Up_h_Slen.Fill(Slen,w3);
	Up_h_max_diff_sollik.Fill(max_diff_sollik,w3);
	Up_h_diffangle_track_shower.Fill(diffangle_track_shower,w3);
	Up_h_myratio50_cascmuon_over_mu.Fill(myratio50_cascmuon_over_mu,w3);
	Up_h_myratio30_cascmuon_over_mu.Fill(myratio30_cascmuon_over_mu,w3);
	Up_h_myratio50_casc_over_mu.Fill(myratio50_casc_over_mu,w3);
	Up_h_myratio30_casc_over_mu.Fill(myratio30_casc_over_mu,w3);
	Up_h_ratio_closehits_cascmuon_over_mu.Fill(ratio_closehits_cascmuon_over_mu,w3);
	Up_h_ratio_closehits_casc_over_mu.Fill(ratio_closehits_casc_over_mu,w3);
	Up_h_redToT_cascmuon_over_mu.Fill(redToT_cascmuon_over_mu,w3);
	Up_h_redToT_casc_over_mu.Fill(redToT_casc_over_mu,w3);
	Up_h_bjorken_y.Fill(bjorken_y,w3);
	Up_h_diffangle.Fill(diffangle,w3);
	Up_h_LogE_mu.Fill(logE_mu,w3);
	Up_h_LogEbundle.Fill(logEbundle,w3);
	//Up_//h_logE_mu_max.Fill(TMath::Log10(E_mu_max),w3);
	Up_h_cos_zen_mu.Fill(cos_zen_mu,w3);
	//------------------------------------------------------
	
        Up_Lik_vs_NNhits.Fill(Nhits,jlik,w3);
	Up_Lik_vs_Q1value.Fill(Q1,jlik,w3);
	Up_NNhits_vs_Q1value.Fill(Q1,Nhits,w3);
	Up_Nlines_vs_zenith.Fill(zenith,num_triggered_lines,w3);
	Up_Chere_Nlines_vs_zenith.Fill(zenith,num_cherenkov_lines,w3);
	Up_TrLengthIT_2_vs_Slen.Fill(Slen,TrLengthIT_2,w3);
	Up_Diffangle_vs_run.Fill(pseudo_runid,diffangle,w3);
	Up_LogEreco_vs_CosZen.Fill(cos_zen,logEreco,w3);
	Up_LogEreco_vs_LogEreco2.Fill(logEreco2,logEreco,w3);
	Up_LogEreco_vs_LogE_mu.Fill(logE_mu,logEreco,w3);
	Up_LogEreco2_vs_LogE_mu.Fill(logE_mu,logEreco2,w3);
	Up_LogEreco_vs_LogE_neu.Fill(logE_nu,logEreco,w3);
	Up_LogEreco2_vs_LogE_neu.Fill(logE_nu,logEreco2,w3);
	Up_LogEresolution_mu_vs_CosZen.Fill(cos_zen,logEresolution_mu,w3);
	Up_LogEresolution_mu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_mu,w3);
	Up_LogEresolution_neu_vs_CosZen.Fill(cos_zen,logEresolution_neu,w3);
	Up_LogEresolution_neu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_neu,w3);


      }//end up



      //Anti-noise cuts
      if(anti_noise_cuts){
	C_CosZen.Fill(cos_zen,w3);
	C_Lik.Fill(jlik,w3);
	C_LogBeta0.Fill(logbeta0,w3);
	C_Beta0_deg.Fill(logbeta0_deg,w3);
	C_DOMs.Fill(num_triggered_doms,w3);
	C_Cherenkov_DOMs.Fill(num_cherenkov_doms,w3);
	C_Q1value.Fill(Q1,w3);
	C_NNhits.Fill(Nhits,w3);
	C_LogEreco.Fill(logEreco,w3);
	C_TrLen.Fill(TrLengthIT_2,w3);
	C_Npe.Fill(SnpeT,w3);
	C_Ratiopmts.Fill(ratio_cherenkov_pmts,w3);
	C_Ratiodoms.Fill(ratio_cherenkov_doms,w3);
	C_Rvalue.Fill(R,w3);
	C_All_triggered_hits.Fill(num_triggered_hits,w3);
	C_Phi.Fill(fi,w3);
	C_Zenith_shower.Fill(zenith_shower,w3);
	C_DLik.Fill(dlik,w3);
	C_Ratio6.Fill(ratio6,w3);
	C_MaxZenSol.Fill(max_zen_sol,w3);
	C_RedToT_muon.Fill(redToT_muon,w3);
	C_Zenith.Fill(zenith,w3);
	C_Delta_ZeniTh_track_shower.Fill(delta_zenith_track_shower,w3);
	C_MinZenSol.Fill(min_zen_sol,w3);
	C_ItoLen.Fill(itoverlen,w3);
	C_Ratio330.Fill(ratio330,w3);
	C_Num_good_sol.Fill(num_of_good_sol,w3);
	C_PosZ.Fill(best_trk_pos_z,w3);
	C_LogEnu.Fill(logE_nu,w3);
	C_LogAEnergy.Fill(logEreco_shower,w3);
	C_NLines.Fill(num_triggered_lines,w3);
	C_LogEreco_cor.Fill(logEreco2,w3);
	C_Chere_NLines.Fill(num_cherenkov_lines,w3);
	C_Ratio_upper_hemisphere.Fill(Nhits_upper/(Nhits_upper + Nhits_lower), w3);
	C_Ratio_lower_hemisphere.Fill(Nhits_lower/(Nhits_upper + Nhits_lower), w3);
	C_Ratio_border_upper_hemisphere.Fill(Nhits_border_upper/Nborder_hits, w3);
	C_Ratio_border_lower_hemisphere.Fill(Nhits_border_lower/Nborder_hits, w3);
	C_Ratio1.Fill(ratio1,w3);
	C_Ratio110.Fill(ratio110,w3);
	C_Ratio130.Fill(ratio130,w3);
	C_NtrackIT_ratio.Fill(NtrackIT/Nhits, w3);
	C_TrLenIT_3.Fill(TrLengthIT_3, w3);

	C_h_ToT_border_mu.Fill(ToT_border_mu,w3);
	C_h_ToT_border_casc.Fill(ToT_border_casc,w3);
	C_h_ToT_border_cascmu.Fill(ToT_border_cascmu,w3);
	C_h_ToT_trig.Fill(ToT_trig,w3);
	C_h_max_ToT_trig.Fill(max_ToT_trig,w3);
	C_h_ToT_IT.Fill(ToT_IT,w3);
	C_h_ToT_allIT.Fill(ToT_allIT,w3);
	C_h_sum_ToT_casc.Fill(sum_ToT_casc,w3);
	C_h_sum_ToT_mu.Fill(sum_ToT_mu,w3);
	C_h_sum_ToT_casc_mu.Fill(sum_ToT_casc_mu,w3);

	//-----------------------------------------------
	C_h_Ratio430.Fill(ratio430,w3);
	C_h_Ratio5.Fill(ratio5,w3);
	C_h_DiffDistance_mu.Fill(diff_dist_mu,w3);
	C_h_RatioCloseHits_mu.Fill(ratio_closehits_muon,w3);
	C_h_NtrackEarly.Fill(NtrackEarly,w3);
	C_h_D.Fill(D,w3);
	C_h_NtrackIT30.Fill(NtrackIT30,w3);
	C_h_NtrackLate.Fill(NtrackLate,w3);
	C_h_myratio50_muon.Fill(myratio50_muon,w3);
	C_h_ratio_cherenkov_lines.Fill(ratio_cherenkov_lines,w3);
	C_h_SnpeT.Fill(SnpeT,w3);
	C_h_mean_tres_it.Fill(mean_tres_it,w3);
	C_h_max_lik_up.Fill(max_lik_up,w3);
	C_h_max_lik_down.Fill(max_lik_down,w3);
	C_h_diff_theta.Fill(diff_theta,w3);
	C_h_diff_dist_casc_mu.Fill(diff_dist_casc_mu,w3);
	C_h_diff_dist_casc.Fill(diff_dist_casc,w3);
	C_h_ratio_closehits_cascmuon.Fill(ratio_closehits_cascmuon,w3);
	C_h_ratio_closehits_casc.Fill(ratio_closehits_casc,w3);
	C_h_redToT_cascmuon.Fill(redToT_cascmuon,w3);
	C_h_redToT_casc.Fill(redToT_casc,w3);
	C_h_myratio50_cascmuon.Fill(myratio50_cascmuon,w3);
	C_h_myratio50_casc.Fill(myratio50_casc,w3);
	C_h_myratio30_cascmuon.Fill(myratio30_cascmuon,w3);
	C_h_myratio30_muon.Fill(myratio30_muon,w3);
	C_h_myratio30_casc.Fill(myratio30_casc,w3);
	C_h_min_diff_sollik.Fill(min_diff_sollik,w3);
	C_h_beta0_shower_deg.Fill(beta0_shower_deg,w3);
	C_h_lik_shower.Fill(lik_shower,w3);
	C_h_best_trk_pos_shower_z.Fill(best_trk_pos_shower_z,w3);
	C_h_normdlik.Fill(normdlik,w3);
	C_h_upsol.Fill(upsol,w3);
	C_h_Slen.Fill(Slen,w3);
	C_h_max_diff_sollik.Fill(max_diff_sollik,w3);
	C_h_diffangle_track_shower.Fill(diffangle_track_shower,w3);
	C_h_myratio50_cascmuon_over_mu.Fill(myratio50_cascmuon_over_mu,w3);
	C_h_myratio30_cascmuon_over_mu.Fill(myratio30_cascmuon_over_mu,w3);
	C_h_myratio50_casc_over_mu.Fill(myratio50_casc_over_mu,w3);
	C_h_myratio30_casc_over_mu.Fill(myratio30_casc_over_mu,w3);
	C_h_ratio_closehits_cascmuon_over_mu.Fill(ratio_closehits_cascmuon_over_mu,w3);
	C_h_ratio_closehits_casc_over_mu.Fill(ratio_closehits_casc_over_mu,w3);
	C_h_redToT_cascmuon_over_mu.Fill(redToT_cascmuon_over_mu,w3);
	C_h_redToT_casc_over_mu.Fill(redToT_casc_over_mu,w3);
	C_h_bjorken_y.Fill(bjorken_y,w3);
	C_h_diffangle.Fill(diffangle,w3);
	C_h_LogE_mu.Fill(logE_mu,w3);
	C_h_LogEbundle.Fill(logEbundle,w3);
	//c_//h_logE_mu_max.Fill(TMath::Log10(E_mu_max),w3);
	C_h_cos_zen_mu.Fill(cos_zen_mu,w3);
	//----------------------------------------------

	C_Lik_vs_NNhits.Fill(Nhits,jlik,w3);
	C_Lik_vs_Q1value.Fill(Q1,jlik,w3);
	C_NNhits_vs_Q1value.Fill(Q1,Nhits,w3);
	C_Nlines_vs_zenith.Fill(zenith,num_triggered_lines,w3);
	C_Chere_Nlines_vs_zenith.Fill(zenith,num_cherenkov_lines,w3);
	C_TrLengthIT_2_vs_Slen.Fill(Slen,TrLengthIT_2,w3);
	C_Diffangle_vs_run.Fill(pseudo_runid,diffangle,w3);
	C_LogEreco_vs_CosZen.Fill(cos_zen,logEreco,w3);	
	C_LogEreco_vs_LogEreco2.Fill(logEreco2,logEreco,w3);
	C_LogEreco_vs_LogE_mu.Fill(logE_mu,logEreco,w3);
	C_LogEreco2_vs_LogE_mu.Fill(logE_mu,logEreco2,w3);
	C_LogEreco_vs_LogE_neu.Fill(logE_nu,logEreco,w3);
	C_LogEreco2_vs_LogE_neu.Fill(logE_nu,logEreco2,w3);
	C_LogEresolution_mu_vs_CosZen.Fill(cos_zen,logEresolution_mu,w3);
	C_LogEresolution_mu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_mu,w3);
	C_LogEresolution_neu_vs_CosZen.Fill(cos_zen,logEresolution_neu,w3);
	C_LogEresolution_neu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_neu,w3);

	//up
	if(up_cut){
	  C_Up_CosZen.Fill(cos_zen,w3);
	  C_Up_Lik.Fill(jlik,w3);
	  C_Up_LogBeta0.Fill(logbeta0,w3);
	  C_Up_Beta0_deg.Fill(logbeta0_deg,w3);
	  C_Up_DOMs.Fill(num_triggered_doms,w3);
	  C_Up_Cherenkov_DOMs.Fill(num_cherenkov_doms,w3);
	  C_Up_Q1value.Fill(Q1,w3);
	  C_Up_NNhits.Fill(Nhits,w3);
	  C_Up_LogEreco.Fill(logEreco,w3);
	  C_Up_TrLen.Fill(TrLengthIT_2,w3);
	  C_Up_Npe.Fill(SnpeT,w3);
	  C_Up_Ratiopmts.Fill(ratio_cherenkov_pmts,w3);
	  C_Up_Ratiodoms.Fill(ratio_cherenkov_doms,w3);
	  C_Up_Rvalue.Fill(R,w3);
	  C_Up_All_triggered_hits.Fill(num_triggered_hits,w3);
	  C_Up_Phi.Fill(fi,w3);
	  C_Up_Zenith_shower.Fill(zenith_shower,w3);
	  C_Up_DLik.Fill(dlik,w3);
	  C_Up_Ratio6.Fill(ratio6,w3);
	  C_Up_MaxZenSol.Fill(max_zen_sol,w3);
	  C_Up_RedToT_muon.Fill(redToT_muon,w3);
	  C_Up_Zenith.Fill(zenith,w3);
	  C_Up_Delta_ZeniTh_track_shower.Fill(delta_zenith_track_shower,w3);
	  C_Up_MinZenSol.Fill(min_zen_sol,w3);
	  C_Up_ItoLen.Fill(itoverlen,w3);
	  C_Up_Ratio330.Fill(ratio330,w3);
	  C_Up_Num_good_sol.Fill(num_of_good_sol,w3);
	  C_Up_PosZ.Fill(best_trk_pos_z,w3);
	  C_Up_LogEnu.Fill(logE_nu,w3);
	  C_Up_LogAEnergy.Fill(logEreco_shower,w3);
	  C_Up_NLines.Fill(num_triggered_lines,w3);
	  C_Up_LogEreco_cor.Fill(logEreco2,w3);
	  C_Up_Chere_NLines.Fill(num_cherenkov_lines,w3);
	  C_Up_Ratio_upper_hemisphere.Fill(Nhits_upper/(Nhits_upper + Nhits_lower), w3);
	  C_Up_Ratio_lower_hemisphere.Fill(Nhits_lower/(Nhits_upper + Nhits_lower), w3);
	  C_Up_Ratio_border_upper_hemisphere.Fill(Nhits_border_upper/Nborder_hits, w3);
	  C_Up_Ratio_border_lower_hemisphere.Fill(Nhits_border_lower/Nborder_hits, w3);
	  C_Up_Ratio1.Fill(ratio1,w3);
	  C_Up_Ratio110.Fill(ratio110,w3);
	  C_Up_Ratio130.Fill(ratio130,w3);
	  C_Up_NtrackIT_ratio.Fill(NtrackIT/Nhits, w3);
	  C_Up_TrLenIT_3.Fill(TrLengthIT_3, w3);

	  C_Up_h_ToT_border_mu.Fill(ToT_border_mu,w3);
	  C_Up_h_ToT_border_casc.Fill(ToT_border_casc,w3);
	  C_Up_h_ToT_border_cascmu.Fill(ToT_border_cascmu,w3);
	  C_Up_h_ToT_trig.Fill(ToT_trig,w3);
	  C_Up_h_max_ToT_trig.Fill(max_ToT_trig,w3);
	  C_Up_h_ToT_IT.Fill(ToT_IT,w3);
	  C_Up_h_ToT_allIT.Fill(ToT_allIT,w3);
	  C_Up_h_sum_ToT_casc.Fill(sum_ToT_casc,w3);
	  C_Up_h_sum_ToT_mu.Fill(sum_ToT_mu,w3);
	  C_Up_h_sum_ToT_casc_mu.Fill(sum_ToT_casc_mu,w3);

	  //-----------------------------------------------
	  C_Up_h_Ratio430.Fill(ratio430,w3);
	  C_Up_h_Ratio5.Fill(ratio5,w3);
	  C_Up_h_DiffDistance_mu.Fill(diff_dist_mu,w3);
	  C_Up_h_RatioCloseHits_mu.Fill(ratio_closehits_muon,w3);
	  C_Up_h_NtrackEarly.Fill(NtrackEarly,w3);
	  C_Up_h_D.Fill(D,w3);
	  C_Up_h_NtrackIT30.Fill(NtrackIT30,w3);
	  C_Up_h_NtrackLate.Fill(NtrackLate,w3);
	  C_Up_h_myratio50_muon.Fill(myratio50_muon,w3);
	  C_Up_h_ratio_cherenkov_lines.Fill(ratio_cherenkov_lines,w3);
	  C_Up_h_SnpeT.Fill(SnpeT,w3);
	  C_Up_h_mean_tres_it.Fill(mean_tres_it,w3);
	  C_Up_h_max_lik_up.Fill(max_lik_up,w3);
	  C_Up_h_max_lik_down.Fill(max_lik_down,w3);
	  C_Up_h_diff_theta.Fill(diff_theta,w3);
	  C_Up_h_diff_dist_casc_mu.Fill(diff_dist_casc_mu,w3);
	  C_Up_h_diff_dist_casc.Fill(diff_dist_casc,w3);
	  C_Up_h_ratio_closehits_cascmuon.Fill(ratio_closehits_cascmuon,w3);
	  C_Up_h_ratio_closehits_casc.Fill(ratio_closehits_casc,w3);
	  C_Up_h_redToT_cascmuon.Fill(redToT_cascmuon,w3);
	  C_Up_h_redToT_casc.Fill(redToT_casc,w3);
	  C_Up_h_myratio50_cascmuon.Fill(myratio50_cascmuon,w3);
	  C_Up_h_myratio50_casc.Fill(myratio50_casc,w3);
	  C_Up_h_myratio30_cascmuon.Fill(myratio30_cascmuon,w3);
	  C_Up_h_myratio30_muon.Fill(myratio30_muon,w3);
	  C_Up_h_myratio30_casc.Fill(myratio30_casc,w3);
	  C_Up_h_min_diff_sollik.Fill(min_diff_sollik,w3);
	  C_Up_h_beta0_shower_deg.Fill(beta0_shower_deg,w3);
	  C_Up_h_lik_shower.Fill(lik_shower,w3);
	  C_Up_h_best_trk_pos_shower_z.Fill(best_trk_pos_shower_z,w3);
	  C_Up_h_normdlik.Fill(normdlik,w3);
	  C_Up_h_upsol.Fill(upsol,w3);
	  C_Up_h_Slen.Fill(Slen,w3);
	  C_Up_h_max_diff_sollik.Fill(max_diff_sollik,w3);
	  C_Up_h_diffangle_track_shower.Fill(diffangle_track_shower,w3);
	  C_Up_h_myratio50_cascmuon_over_mu.Fill(myratio50_cascmuon_over_mu,w3);
	  C_Up_h_myratio30_cascmuon_over_mu.Fill(myratio30_cascmuon_over_mu,w3);
	  C_Up_h_myratio50_casc_over_mu.Fill(myratio50_casc_over_mu,w3);
	  C_Up_h_myratio30_casc_over_mu.Fill(myratio30_casc_over_mu,w3);
	  C_Up_h_ratio_closehits_cascmuon_over_mu.Fill(ratio_closehits_cascmuon_over_mu,w3);
	  C_Up_h_ratio_closehits_casc_over_mu.Fill(ratio_closehits_casc_over_mu,w3);
	  C_Up_h_redToT_cascmuon_over_mu.Fill(redToT_cascmuon_over_mu,w3);
	  C_Up_h_redToT_casc_over_mu.Fill(redToT_casc_over_mu,w3);
	  C_Up_h_bjorken_y.Fill(bjorken_y,w3);
	  C_Up_h_diffangle.Fill(diffangle,w3);
	  C_Up_h_LogE_mu.Fill(logE_mu,w3);
	  C_Up_h_LogEbundle.Fill(logEbundle,w3);
	  //C_Up_//h_logE_mu_max.Fill(TMath::Log10(E_mu_max),w3);
	  C_Up_h_cos_zen_mu.Fill(cos_zen_mu,w3);
	  //----------------------------------------------


	  C_Up_Lik_vs_NNhits.Fill(Nhits,jlik,w3);
	  C_Up_Lik_vs_Q1value.Fill(Q1,jlik,w3);
	  C_Up_NNhits_vs_Q1value.Fill(Q1,Nhits,w3);
	  C_Up_Nlines_vs_zenith.Fill(zenith,num_triggered_lines,w3);
	  C_Up_Chere_Nlines_vs_zenith.Fill(zenith,num_cherenkov_lines,w3);
	  C_Up_TrLengthIT_2_vs_Slen.Fill(Slen,TrLengthIT_2,w3);
	  C_Up_Diffangle_vs_run.Fill(pseudo_runid,diffangle,w3);
	  C_Up_LogEreco_vs_CosZen.Fill(cos_zen,logEreco,w3);
	  C_Up_LogEreco_vs_LogEreco2.Fill(logEreco2,logEreco,w3);
	  C_Up_LogEreco_vs_LogE_mu.Fill(logE_mu,logEreco,w3);
	  C_Up_LogEreco2_vs_LogE_mu.Fill(logE_mu,logEreco2,w3);
	  C_Up_LogEreco_vs_LogE_neu.Fill(logE_nu,logEreco,w3);
	  C_Up_LogEreco2_vs_LogE_neu.Fill(logE_nu,logEreco2,w3);
	  C_Up_LogEresolution_mu_vs_CosZen.Fill(cos_zen,logEresolution_mu,w3);
	  C_Up_LogEresolution_mu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_mu,w3);
	  C_Up_LogEresolution_neu_vs_CosZen.Fill(cos_zen,logEresolution_neu,w3);
	  C_Up_LogEresolution_neu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_neu,w3);


	}//end up
      }//end anti-noise

      //Frst cuts
      if(first_precuts){
	First_CosZen.Fill(cos_zen,w3);
	First_Lik.Fill(jlik,w3);
	First_LogBeta0.Fill(logbeta0,w3);
	First_Beta0_deg.Fill(logbeta0_deg,w3);
	First_DOMs.Fill(num_triggered_doms,w3);
	First_Cherenkov_DOMs.Fill(num_cherenkov_doms,w3);
	First_Q1value.Fill(Q1,w3);
	First_NNhits.Fill(Nhits,w3);
	First_LogEreco.Fill(logEreco,w3);    
	First_TrLen.Fill(TrLengthIT_2,w3);
	First_Npe.Fill(SnpeT,w3);
	First_Ratiopmts.Fill(ratio_cherenkov_pmts,w3);
	First_Ratiodoms.Fill(ratio_cherenkov_doms,w3);
	First_Rvalue.Fill(R,w3);
	First_All_triggered_hits.Fill(num_triggered_hits,w3);
	First_Phi.Fill(fi,w3);
	First_Zenith_shower.Fill(zenith_shower,w3);
	First_DLik.Fill(dlik,w3);
	First_Ratio6.Fill(ratio6,w3);
	First_MaxZenSol.Fill(max_zen_sol,w3);
	First_RedToT_muon.Fill(redToT_muon,w3);
	First_Zenith.Fill(zenith,w3);
	First_Delta_ZeniTh_track_shower.Fill(delta_zenith_track_shower,w3);
	First_MinZenSol.Fill(min_zen_sol,w3);
	First_ItoLen.Fill(itoverlen,w3);
	First_Ratio330.Fill(ratio330,w3);
	First_Num_good_sol.Fill(num_of_good_sol,w3);
	First_PosZ.Fill(best_trk_pos_z,w3);
	First_LogEnu.Fill(logE_nu,w3);
	First_LogAEnergy.Fill(logEreco_shower,w3);
	First_NLines.Fill(num_triggered_lines,w3);
	First_LogEreco_cor.Fill(logEreco2,w3);
	First_Chere_NLines.Fill(num_cherenkov_lines,w3);
	First_Ratio_upper_hemisphere.Fill(Nhits_upper/(Nhits_upper + Nhits_lower), w3);
	First_Ratio_lower_hemisphere.Fill(Nhits_lower/(Nhits_upper + Nhits_lower), w3);
	First_Ratio_border_upper_hemisphere.Fill(Nhits_border_upper/Nborder_hits, w3);
	First_Ratio_border_lower_hemisphere.Fill(Nhits_border_lower/Nborder_hits, w3);
	First_Ratio1.Fill(ratio1,w3);
	First_Ratio110.Fill(ratio110,w3);
	First_Ratio130.Fill(ratio130,w3);
	First_NtrackIT_ratio.Fill(NtrackIT/Nhits, w3);
	First_TrLenIT_3.Fill(TrLengthIT_3, w3);

	First_h_ToT_border_mu.Fill(ToT_border_mu,w3);
	First_h_ToT_border_casc.Fill(ToT_border_casc,w3);
	First_h_ToT_border_cascmu.Fill(ToT_border_cascmu,w3);
	First_h_ToT_trig.Fill(ToT_trig,w3);
	First_h_max_ToT_trig.Fill(max_ToT_trig,w3);
	First_h_ToT_IT.Fill(ToT_IT,w3);
	First_h_ToT_allIT.Fill(ToT_allIT,w3);
	First_h_sum_ToT_casc.Fill(sum_ToT_casc,w3);
	First_h_sum_ToT_mu.Fill(sum_ToT_mu,w3);
	First_h_sum_ToT_casc_mu.Fill(sum_ToT_casc_mu,w3);

	//-----------------------------------------------
	First_h_Ratio430.Fill(ratio430,w3);
	First_h_Ratio5.Fill(ratio5,w3);
	First_h_DiffDistance_mu.Fill(diff_dist_mu,w3);
	First_h_RatioCloseHits_mu.Fill(ratio_closehits_muon,w3);
	First_h_NtrackEarly.Fill(NtrackEarly,w3);
	First_h_D.Fill(D,w3);
	First_h_NtrackIT30.Fill(NtrackIT30,w3);
	First_h_NtrackLate.Fill(NtrackLate,w3);
	First_h_myratio50_muon.Fill(myratio50_muon,w3);
	First_h_ratio_cherenkov_lines.Fill(ratio_cherenkov_lines,w3);
	First_h_SnpeT.Fill(SnpeT,w3);
	First_h_mean_tres_it.Fill(mean_tres_it,w3);
	First_h_max_lik_up.Fill(max_lik_up,w3);
	First_h_max_lik_down.Fill(max_lik_down,w3);
	First_h_diff_theta.Fill(diff_theta,w3);
	First_h_diff_dist_casc_mu.Fill(diff_dist_casc_mu,w3);
	First_h_diff_dist_casc.Fill(diff_dist_casc,w3);
	First_h_ratio_closehits_cascmuon.Fill(ratio_closehits_cascmuon,w3);
	First_h_ratio_closehits_casc.Fill(ratio_closehits_casc,w3);
	First_h_redToT_cascmuon.Fill(redToT_cascmuon,w3);
	First_h_redToT_casc.Fill(redToT_casc,w3);
	First_h_myratio50_cascmuon.Fill(myratio50_cascmuon,w3);
	First_h_myratio50_casc.Fill(myratio50_casc,w3);
	First_h_myratio30_cascmuon.Fill(myratio30_cascmuon,w3);
	First_h_myratio30_muon.Fill(myratio30_muon,w3);
	First_h_myratio30_casc.Fill(myratio30_casc,w3);
	First_h_min_diff_sollik.Fill(min_diff_sollik,w3);
	First_h_beta0_shower_deg.Fill(beta0_shower_deg,w3);
	First_h_lik_shower.Fill(lik_shower,w3);
	First_h_best_trk_pos_shower_z.Fill(best_trk_pos_shower_z,w3);
	First_h_normdlik.Fill(normdlik,w3);
	First_h_upsol.Fill(upsol,w3);
	First_h_Slen.Fill(Slen,w3);
	First_h_max_diff_sollik.Fill(max_diff_sollik,w3);
	First_h_diffangle_track_shower.Fill(diffangle_track_shower,w3);
	First_h_myratio50_cascmuon_over_mu.Fill(myratio50_cascmuon_over_mu,w3);
	First_h_myratio30_cascmuon_over_mu.Fill(myratio30_cascmuon_over_mu,w3);
	First_h_myratio50_casc_over_mu.Fill(myratio50_casc_over_mu,w3);
	First_h_myratio30_casc_over_mu.Fill(myratio30_casc_over_mu,w3);
	First_h_ratio_closehits_cascmuon_over_mu.Fill(ratio_closehits_cascmuon_over_mu,w3);
	First_h_ratio_closehits_casc_over_mu.Fill(ratio_closehits_casc_over_mu,w3);
	First_h_redToT_cascmuon_over_mu.Fill(redToT_cascmuon_over_mu,w3);
	First_h_redToT_casc_over_mu.Fill(redToT_casc_over_mu,w3);
	First_h_bjorken_y.Fill(bjorken_y,w3);
	First_h_diffangle.Fill(diffangle,w3);
	First_h_LogE_mu.Fill(logE_mu,w3);
	First_h_LogEbundle.Fill(logEbundle,w3);
	//First_//h_logE_mu_max.Fill(TMath::Log10(E_mu_max),w3);
	First_h_cos_zen_mu.Fill(cos_zen_mu,w3);
	//----------------------------------------------
      
	First_Lik_vs_NNhits.Fill(Nhits,jlik,w3);
	First_Lik_vs_Q1value.Fill(Q1,jlik,w3);
	First_NNhits_vs_Q1value.Fill(Q1,Nhits,w3);
	First_Nlines_vs_zenith.Fill(zenith,num_triggered_lines,w3);
	First_Chere_Nlines_vs_zenith.Fill(zenith,num_cherenkov_lines,w3);
	First_TrLengthIT_2_vs_Slen.Fill(Slen,TrLengthIT_2,w3);
	First_Diffangle_vs_run.Fill(pseudo_runid,diffangle,w3);
	First_LogEreco_vs_CosZen.Fill(cos_zen,logEreco,w3);
	First_LogEreco_vs_LogEreco2.Fill(logEreco2,logEreco,w3);
	First_LogEreco_vs_LogE_mu.Fill(logE_mu,logEreco,w3);
	First_LogEreco2_vs_LogE_mu.Fill(logE_mu,logEreco2,w3);
	First_LogEreco_vs_LogE_neu.Fill(logE_nu,logEreco,w3);
	First_LogEreco2_vs_LogE_neu.Fill(logE_nu,logEreco2,w3);
	First_LogEresolution_mu_vs_CosZen.Fill(cos_zen,logEresolution_mu,w3);
	First_LogEresolution_mu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_mu,w3);
	First_LogEresolution_neu_vs_CosZen.Fill(cos_zen,logEresolution_neu,w3);
	First_LogEresolution_neu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_neu,w3);


	//up
	if(up_cut){
	  First_Up_CosZen.Fill(cos_zen,w3);
	  First_Up_Lik.Fill(jlik,w3);
	  First_Up_LogBeta0.Fill(logbeta0,w3);
	  First_Up_Beta0_deg.Fill(logbeta0_deg,w3);
	  First_Up_DOMs.Fill(num_triggered_doms,w3);
	  First_Up_Cherenkov_DOMs.Fill(num_cherenkov_doms,w3);
	  First_Up_Q1value.Fill(Q1,w3);
	  First_Up_NNhits.Fill(Nhits,w3);
	  First_Up_LogEreco.Fill(logEreco,w3);
	  First_Up_TrLen.Fill(TrLengthIT_2,w3);
	  First_Up_Npe.Fill(SnpeT,w3);
	  First_Up_Ratiopmts.Fill(ratio_cherenkov_pmts,w3);
	  First_Up_Ratiodoms.Fill(ratio_cherenkov_doms,w3);
	  First_Up_Rvalue.Fill(R,w3);
	  First_Up_All_triggered_hits.Fill(num_triggered_hits,w3);
	  First_Up_Phi.Fill(fi,w3);
	  First_Up_Zenith_shower.Fill(zenith_shower,w3);
	  First_Up_DLik.Fill(dlik,w3);
	  First_Up_Ratio6.Fill(ratio6,w3);
	  First_Up_MaxZenSol.Fill(max_zen_sol,w3);
	  First_Up_RedToT_muon.Fill(redToT_muon,w3);
	  First_Up_Zenith.Fill(zenith,w3);
	  First_Up_Delta_ZeniTh_track_shower.Fill(delta_zenith_track_shower,w3);
	  First_Up_MinZenSol.Fill(min_zen_sol,w3);
	  First_Up_ItoLen.Fill(itoverlen,w3);
	  First_Up_Ratio330.Fill(ratio330,w3);
	  First_Up_Num_good_sol.Fill(num_of_good_sol,w3);
	  First_Up_PosZ.Fill(best_trk_pos_z,w3);
	  First_Up_LogEnu.Fill(logE_nu,w3);
	  First_Up_LogAEnergy.Fill(logEreco_shower,w3);
	  First_Up_NLines.Fill(num_triggered_lines,w3);
	  First_Up_LogEreco_cor.Fill(logEreco2,w3);
	  First_Up_Chere_NLines.Fill(num_cherenkov_lines,w3);
	  First_Up_Ratio_upper_hemisphere.Fill(Nhits_upper/(Nhits_upper + Nhits_lower), w3);
	  First_Up_Ratio_lower_hemisphere.Fill(Nhits_lower/(Nhits_upper + Nhits_lower), w3);
	  First_Up_Ratio_border_upper_hemisphere.Fill(Nhits_border_upper/Nborder_hits, w3);
	  First_Up_Ratio_border_lower_hemisphere.Fill(Nhits_border_lower/Nborder_hits, w3);
	  First_Up_Ratio1.Fill(ratio1,w3);
	  First_Up_Ratio110.Fill(ratio110,w3);
	  First_Up_Ratio130.Fill(ratio130,w3);
	  First_Up_NtrackIT_ratio.Fill(NtrackIT/Nhits, w3);
	  First_Up_TrLenIT_3.Fill(TrLengthIT_3, w3);

	  First_Up_h_ToT_border_mu.Fill(ToT_border_mu,w3);
	  First_Up_h_ToT_border_casc.Fill(ToT_border_casc,w3);
	  First_Up_h_ToT_border_cascmu.Fill(ToT_border_cascmu,w3);
	  First_Up_h_ToT_trig.Fill(ToT_trig,w3);
	  First_Up_h_max_ToT_trig.Fill(max_ToT_trig,w3);
	  First_Up_h_ToT_IT.Fill(ToT_IT,w3);
	  First_Up_h_ToT_allIT.Fill(ToT_allIT,w3);
	  First_Up_h_sum_ToT_casc.Fill(sum_ToT_casc,w3);
	  First_Up_h_sum_ToT_mu.Fill(sum_ToT_mu,w3);
	  First_Up_h_sum_ToT_casc_mu.Fill(sum_ToT_casc_mu,w3);

	  //-----------------------------------------------
	  First_Up_h_Ratio430.Fill(ratio430,w3);
	  First_Up_h_Ratio5.Fill(ratio5,w3);
	  First_Up_h_DiffDistance_mu.Fill(diff_dist_mu,w3);
	  First_Up_h_RatioCloseHits_mu.Fill(ratio_closehits_muon,w3);
	  First_Up_h_NtrackEarly.Fill(NtrackEarly,w3);
	  First_Up_h_D.Fill(D,w3);
	  First_Up_h_NtrackIT30.Fill(NtrackIT30,w3);
	  First_Up_h_NtrackLate.Fill(NtrackLate,w3);
	  First_Up_h_myratio50_muon.Fill(myratio50_muon,w3);
	  First_Up_h_ratio_cherenkov_lines.Fill(ratio_cherenkov_lines,w3);
	  First_Up_h_SnpeT.Fill(SnpeT,w3);
	  First_Up_h_mean_tres_it.Fill(mean_tres_it,w3);
	  First_Up_h_max_lik_up.Fill(max_lik_up,w3);
	  First_Up_h_max_lik_down.Fill(max_lik_down,w3);
	  First_Up_h_diff_theta.Fill(diff_theta,w3);
	  First_Up_h_diff_dist_casc_mu.Fill(diff_dist_casc_mu,w3);
	  First_Up_h_diff_dist_casc.Fill(diff_dist_casc,w3);
	  First_Up_h_ratio_closehits_cascmuon.Fill(ratio_closehits_cascmuon,w3);
	  First_Up_h_ratio_closehits_casc.Fill(ratio_closehits_casc,w3);
	  First_Up_h_redToT_cascmuon.Fill(redToT_cascmuon,w3);
	  First_Up_h_redToT_casc.Fill(redToT_casc,w3);
	  First_Up_h_myratio50_cascmuon.Fill(myratio50_cascmuon,w3);
	  First_Up_h_myratio50_casc.Fill(myratio50_casc,w3);
	  First_Up_h_myratio30_cascmuon.Fill(myratio30_cascmuon,w3);
	  First_Up_h_myratio30_muon.Fill(myratio30_muon,w3);
	  First_Up_h_myratio30_casc.Fill(myratio30_casc,w3);
	  First_Up_h_min_diff_sollik.Fill(min_diff_sollik,w3);
	  First_Up_h_beta0_shower_deg.Fill(beta0_shower_deg,w3);
	  First_Up_h_lik_shower.Fill(lik_shower,w3);
	  First_Up_h_best_trk_pos_shower_z.Fill(best_trk_pos_shower_z,w3);
	  First_Up_h_normdlik.Fill(normdlik,w3);
	  First_Up_h_upsol.Fill(upsol,w3);
	  First_Up_h_Slen.Fill(Slen,w3);
	  First_Up_h_max_diff_sollik.Fill(max_diff_sollik,w3);
	  First_Up_h_diffangle_track_shower.Fill(diffangle_track_shower,w3);
	  First_Up_h_myratio50_cascmuon_over_mu.Fill(myratio50_cascmuon_over_mu,w3);
	  First_Up_h_myratio30_cascmuon_over_mu.Fill(myratio30_cascmuon_over_mu,w3);
	  First_Up_h_myratio50_casc_over_mu.Fill(myratio50_casc_over_mu,w3);
	  First_Up_h_myratio30_casc_over_mu.Fill(myratio30_casc_over_mu,w3);
	  First_Up_h_ratio_closehits_cascmuon_over_mu.Fill(ratio_closehits_cascmuon_over_mu,w3);
	  First_Up_h_ratio_closehits_casc_over_mu.Fill(ratio_closehits_casc_over_mu,w3);
	  First_Up_h_redToT_cascmuon_over_mu.Fill(redToT_cascmuon_over_mu,w3);
	  First_Up_h_redToT_casc_over_mu.Fill(redToT_casc_over_mu,w3);
	  First_Up_h_bjorken_y.Fill(bjorken_y,w3);
	  First_Up_h_diffangle.Fill(diffangle,w3);
	  First_Up_h_LogE_mu.Fill(logE_mu,w3);
	  First_Up_h_LogEbundle.Fill(logEbundle,w3);
	  //First_Up_//h_logE_mu_max.Fill(TMath::Log10(E_mu_max),w3);
	  First_Up_h_cos_zen_mu.Fill(cos_zen_mu,w3);
	  //----------------------------------------------


	  First_Up_Lik_vs_NNhits.Fill(Nhits,jlik,w3);
	  First_Up_Lik_vs_Q1value.Fill(Q1,jlik,w3);
	  First_Up_NNhits_vs_Q1value.Fill(Q1,Nhits,w3);
	  First_Up_Nlines_vs_zenith.Fill(zenith,num_triggered_lines,w3);
	  First_Up_Chere_Nlines_vs_zenith.Fill(zenith,num_cherenkov_lines,w3);
	  First_Up_TrLengthIT_2_vs_Slen.Fill(Slen,TrLengthIT_2,w3);
	  First_Up_Diffangle_vs_run.Fill(pseudo_runid,diffangle,w3);
	  First_Up_LogEreco_vs_CosZen.Fill(cos_zen,logEreco,w3);
	  First_Up_LogEreco_vs_LogEreco2.Fill(logEreco2,logEreco,w3);
	  First_Up_LogEreco_vs_LogE_mu.Fill(logE_mu,logEreco,w3);
	  First_Up_LogEreco2_vs_LogE_mu.Fill(logE_mu,logEreco2,w3);
	  First_Up_LogEreco_vs_LogE_neu.Fill(logE_nu,logEreco,w3);
	  First_Up_LogEreco2_vs_LogE_neu.Fill(logE_nu,logEreco2,w3);
	  First_Up_LogEresolution_mu_vs_CosZen.Fill(cos_zen,logEresolution_mu,w3);
	  First_Up_LogEresolution_mu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_mu,w3);
	  First_Up_LogEresolution_neu_vs_CosZen.Fill(cos_zen,logEresolution_neu,w3);
	  First_Up_LogEresolution_neu_cor_vs_CosZen.Fill(cos_zen,logEresolution_cor_neu,w3);


	}//end up
      }//end first

      //Common Level cut
      if(common_level){
	Common_CosZen.Fill(cos_zen,w3);     
	Common_Lik.Fill(jlik,w3);
	Common_LogBeta0.Fill(logbeta0,w3);
	Common_LogEreco.Fill(logEreco,w3);
	Common_TrLen.Fill(TrLengthIT_2,w3);
	Common_LogEresolution_mu_vs_CosZen.Fill(cos_zen,logEresolution_mu,w3);
	Common_LogEresolution_neu_vs_CosZen.Fill(cos_zen,logEresolution_neu,w3);
      }//end of common level cuts

    }//extra condition

  } //end of loop over entries (over events, that is)
  //======================================================================

  //loop over all runs (QualityParameters tree) NOT events
  double total_livetime=0;
  double run_livetime=0;
  int check_next_run_number=0;
  
   
  for( int k=0 ; k<InfoTree->GetEntries() ; k++ ){
    
    {InfoTree->GetEntry(k+1); check_next_run_number=run_number;} //fetch next "file" run number

    InfoTree->GetEntry(k);

    // ***IMPORTANT NOTE***
    /* If NOT re-declared the values of the TTree's quantitities are  those of the LATEST InfoTree->GetEntry
       e.g. InfoTree->GetEntry(k);
       {InfoTree->GetEntry(k); acurrent_run_number=run_number;}
       {InfoTree->GetEntry(7); check_next_run_number=run_number;}
       livetime will ALWAYS be that of entry 7, EXCEPT if declared before. acurrent_run_number will change normally because it is re-declared each time.
    */
    
    run_livetime += livetime;

    if( run_number != check_next_run_number || k==InfoTree->GetEntries()-1 ){
      
      Allreco_rate.Fill(float(run_number)+0.5,float(reco_events_to_run[run_number])/run_livetime,1.0);
      C_Allreco_rate.Fill(float(run_number)+0.5,float(c_reco_events_to_run[run_number])/run_livetime,1.0);
      First_Allreco_rate.Fill(float(run_number)+0.5,float(first_reco_events_to_run[run_number])/run_livetime,1.0);

      Up_Allreco_rate.Fill(float(run_number)+0.5,float(up_reco_events_to_run[run_number])/run_livetime,1.0);
      C_Up_Allreco_rate.Fill(float(run_number)+0.5,float(c_up_reco_events_to_run[run_number])/run_livetime,1.0);
      First_Up_Allreco_rate.Fill(float(run_number)+0.5,float(first_up_reco_events_to_run[run_number])/run_livetime,1.0);

      Rate3Dmuon.Fill(float(run_number)+0.5,float(evts_3d_muon_to_run[run_number])/run_livetime,1.0);
      Rate3Dshower.Fill(float(run_number)+0.5,float(evts_3d_shower_to_run[run_number])/run_livetime,1.0);
      RateMXshower.Fill(float(run_number)+0.5,float(evts_mx_shower_to_run[run_number])/run_livetime,1.0);
      C_Rate3Dmuon.Fill(float(run_number)+0.5,float(c_evts_3d_muon_to_run[run_number])/run_livetime,1.0);
      C_Rate3Dshower.Fill(float(run_number)+0.5,float(c_evts_3d_shower_to_run[run_number])/run_livetime,1.0);
      C_RateMXshower.Fill(float(run_number)+0.5,float(c_evts_mx_shower_to_run[run_number])/run_livetime,1.0);
      First_Rate3Dmuon.Fill(float(run_number)+0.5,float(first_evts_3d_muon_to_run[run_number])/run_livetime,1.0);
      First_Rate3Dshower.Fill(float(run_number)+0.5,float(first_evts_3d_shower_to_run[run_number])/run_livetime,1.0);
      First_RateMXshower.Fill(float(run_number)+0.5,float(first_evts_mx_shower_to_run[run_number])/run_livetime,1.0);
      Up_Rate3Dmuon.Fill(float(run_number)+0.5,float(up_evts_3d_muon_to_run[run_number])/run_livetime,1.0);
      Up_Rate3Dshower.Fill(float(run_number)+0.5,float(up_evts_3d_shower_to_run[run_number])/run_livetime,1.0);
      Up_RateMXshower.Fill(float(run_number)+0.5,float(up_evts_mx_shower_to_run[run_number])/run_livetime,1.0);
      C_Up_Rate3Dmuon.Fill(float(run_number)+0.5,float(c_up_evts_3d_muon_to_run[run_number])/run_livetime,1.0);
      C_Up_Rate3Dshower.Fill(float(run_number)+0.5,float(c_up_evts_3d_shower_to_run[run_number])/run_livetime,1.0);
      C_Up_RateMXshower.Fill(float(run_number)+0.5,float(c_up_evts_mx_shower_to_run[run_number])/run_livetime,1.0);
      First_Up_Rate3Dmuon.Fill(float(run_number)+0.5,float(first_up_evts_3d_muon_to_run[run_number])/run_livetime,1.0);
      First_Up_Rate3Dshower.Fill(float(run_number)+0.5,float(first_up_evts_3d_shower_to_run[run_number])/run_livetime,1.0);
      First_Up_RateMXshower.Fill(float(run_number)+0.5,float(first_up_evts_mx_shower_to_run[run_number])/run_livetime,1.0);

      std::cout<<"Run: "<<run_number<<" evts: "<<reco_events_to_run[run_number]<<" C_evts: "<<c_reco_events_to_run[run_number]<<" Up_evts: "<<up_reco_events_to_run[run_number]<<" Up_First_evts: "<<first_up_reco_events_to_run[run_number]<<" run_livetime: "<<run_livetime<<std::endl;
      //std::cout<<"Run: "<<run_number<<" run_livetime: "<<run_livetime<<std::endl;
     
      run_livetime=0;

    }//check of current run_number

    total_livetime = total_livetime + livetime;
    
  }//end of loop over runs
  
  std::cout<<" "<<std::endl;
  std::cout<<"-------"<<std::endl;
  std::cout<<"total period's livetime: "<<std::fixed<<std::setprecision(4)<<total_livetime<<" sec"<<std::endl;

  std::cout<<"all events: "<<all_evts<<std::endl;
  std::cout<<"passed events: "<<pass_evts<<std::endl;
  std::cout<<"upgoing passed events: "<<pass_up_evts<<std::endl;
  std::cout<<"anti-noise passed events: "<<pass_C_evts<<std::endl;
  std::cout<<"anti-noise upgoing passed events: "<<pass_C_up_evts<<std::endl;
  std::cout<<"First cuts passed events: "<<pass_First_evts<<std::endl;
  std::cout<<"First cuts upgoing passed events: "<<pass_First_up_evts<<std::endl;

  
  //===== write the analysis histograms ========================================================================
  // OUTPUT
  TFile* f_out = new TFile( RootOut,"RECREATE");

  std::cout<<"I'm writing the output.."<<std::endl;

  //---write 'em-----
  CosZen.Write();//cos_zen
  Lik.Write();//jlik
  LogBeta0.Write();//logbeta0
  Beta0_deg.Write();//logbeta0_deg
  DOMs.Write();//num_triggered_doms
  Cherenkov_DOMs.Write();//num_cherenkov_doms
  Q1value.Write();//Q1
  NNhits.Write();//Nhits
  LogEreco.Write();//logEreco    
  TrLen.Write();//TrLengthIT_2
  Npe.Write();//SnpeT
  Ratiopmts.Write();//ratio_cherenkov_pmts
  Ratiodoms.Write();//ratio_cherenkov_doms
  Rvalue.Write();//R
  All_triggered_hits.Write();//num_triggered_hits
  Phi.Write();//fi
  Zenith_shower.Write();//zenith_shower
  DLik.Write();//dlik
  Ratio6.Write();//ratio6
  MaxZenSol.Write();//max_zen_sol
  RedToT_muon.Write();//redToT_muon
  Zenith.Write();//zenith
  Delta_ZeniTh_track_shower.Write();//delta_zenith_track_shower
  MinZenSol.Write();//min_zen_sol
  ItoLen.Write();//itoverlen
  Ratio330.Write();//ratio330
  Num_good_sol.Write();//num_of_good_sol
  PosZ.Write();//best_trk_pos_z
  LogEnu.Write();//logE
  LogAEnergy.Write();//logeReco_shower
  NLines.Write();//nlines
  LogEreco_cor.Write();
  Chere_NLines.Write();//chere_nlines
  Ratio_upper_hemisphere.Write();
  Ratio_lower_hemisphere.Write();
  Ratio_border_upper_hemisphere.Write();
  Ratio_border_lower_hemisphere.Write();
  Ratio1.Write();
  Ratio110.Write();
  Ratio130.Write();
  NtrackIT_ratio.Write();
  TrLenIT_3.Write();

  h_ToT_border_mu.Write();
  h_ToT_border_casc.Write();
  h_ToT_border_cascmu.Write();
  h_ToT_trig.Write();
  h_max_ToT_trig.Write();
  h_ToT_IT.Write();
  h_ToT_allIT.Write();
  h_sum_ToT_casc.Write();
  h_sum_ToT_mu.Write();
  h_sum_ToT_casc_mu.Write();
  //---------------------------------------------------------------
  h_Ratio430.Write();
  h_Ratio5.Write();
  h_DiffDistance_mu.Write();
  h_RatioCloseHits_mu.Write();
  h_NtrackEarly.Write();
  h_D.Write();
  h_NtrackIT30.Write();
  h_NtrackLate.Write();
  h_myratio50_muon.Write();
  h_ratio_cherenkov_lines.Write();
  h_SnpeT.Write();
  h_mean_tres_it.Write();
  h_max_lik_up.Write();
  h_max_lik_down.Write();
  h_diff_theta.Write();
  h_diff_dist_casc_mu.Write();
  h_diff_dist_casc.Write();
  h_ratio_closehits_cascmuon.Write();
  h_ratio_closehits_casc.Write();
  h_redToT_cascmuon.Write();
  h_redToT_casc.Write();
  h_myratio50_cascmuon.Write();
  h_myratio50_casc.Write();
  h_myratio30_cascmuon.Write();
  h_myratio30_muon.Write();
  h_myratio30_casc.Write();
  h_min_diff_sollik.Write();
  h_beta0_shower_deg.Write();
  h_lik_shower.Write();
  h_best_trk_pos_shower_z.Write();
  h_normdlik.Write();
  h_upsol.Write();
  h_Slen.Write();
  h_max_diff_sollik.Write();
  h_diffangle_track_shower.Write();
  h_myratio50_cascmuon_over_mu.Write();
  h_myratio30_cascmuon_over_mu.Write();
  h_myratio50_casc_over_mu.Write();
  h_myratio30_casc_over_mu.Write();
  h_ratio_closehits_cascmuon_over_mu.Write();
  h_ratio_closehits_casc_over_mu.Write();
  h_redToT_cascmuon_over_mu.Write();
  h_redToT_casc_over_mu.Write();
  h_bjorken_y.Write();
  h_diffangle.Write();
  h_LogE_mu.Write();
  h_LogEbundle.Write();
  //h_logE_mu_max.Write();
  h_cos_zen_mu.Write();
  //-------------------------------------------------------------
      
  Lik_vs_NNhits.Write();//Nhits,jlik
  Lik_vs_Q1value.Write();//Q1,jlik
  NNhits_vs_Q1value.Write();//Q1,Nhits
  Nlines_vs_zenith.Write();//zenith,nlines
  Chere_Nlines_vs_zenith.Write();//zenith,chere_nlines
  TrLengthIT_2_vs_Slen.Write();
  Diffangle_vs_run.Write(); //run vs diffangle (at which run we have each diffangle)
  LogEreco_vs_CosZen.Write();
  LogEreco_vs_LogEreco2.Write();
  LogEreco_vs_LogE_mu.Write();
  LogEreco2_vs_LogE_mu.Write();
  LogEreco_vs_LogE_neu.Write();
  LogEreco2_vs_LogE_neu.Write();
  LogEresolution_mu_vs_CosZen.Write();
  LogEresolution_mu_cor_vs_CosZen.Write();
  LogEresolution_neu_vs_CosZen.Write();
  LogEresolution_neu_cor_vs_CosZen.Write();
     
  Up_CosZen.Write();//cos_zen
  Up_Lik.Write();//jlik
  Up_LogBeta0.Write();//logbeta0
  Up_Beta0_deg.Write();//logbeta0_deg
  Up_DOMs.Write();//num_triggered_doms
  Up_Cherenkov_DOMs.Write();//num_cherenkov_doms
  Up_Q1value.Write();//Q1
  Up_NNhits.Write();//Nhits
  Up_LogEreco.Write();//log10(Ereco)
  Up_TrLen.Write();//TrLengthIT_2
  Up_Npe.Write();//SnpeT
  Up_Ratiopmts.Write();//ratio_cherenkov_pmts
  Up_Ratiodoms.Write();//ratio_cherenkov_doms
  Up_Rvalue.Write();//R
  Up_All_triggered_hits.Write();//num_triggered_hits
  Up_Phi.Write();//fi
  Up_Zenith_shower.Write();//zenith_shower
  Up_DLik.Write();//dlik
  Up_Ratio6.Write();//ratio6
  Up_MaxZenSol.Write();//max_zen_sol
  Up_RedToT_muon.Write();//redToT_muon
  Up_Zenith.Write();//zenith
  Up_Delta_ZeniTh_track_shower.Write();//delta_zenith_track_shower
  Up_MinZenSol.Write();//min_zen_sol
  Up_ItoLen.Write();//itoverlen
  Up_Ratio330.Write();//ratio330
  Up_Num_good_sol.Write();//num_of_good_sol
  Up_PosZ.Write();//best_trk_pos_z
  Up_LogEnu.Write();//logE
  Up_LogAEnergy.Write();//logeReco_shower
  Up_NLines.Write();//nlines
  Up_LogEreco_cor.Write();
  Up_Chere_NLines.Write();//chere_nlines
  Up_Ratio_upper_hemisphere.Write();
  Up_Ratio_lower_hemisphere.Write();
  Up_Ratio_border_upper_hemisphere.Write();
  Up_Ratio_border_lower_hemisphere.Write();
  Up_Ratio1.Write();
  Up_Ratio110.Write();
  Up_Ratio130.Write();
  Up_NtrackIT_ratio.Write();
  Up_TrLenIT_3.Write();

  Up_h_ToT_border_mu.Write();
  Up_h_ToT_border_casc.Write();
  Up_h_ToT_border_cascmu.Write();
  Up_h_ToT_trig.Write();
  Up_h_max_ToT_trig.Write();
  Up_h_ToT_IT.Write();
  Up_h_ToT_allIT.Write();
  Up_h_sum_ToT_casc.Write();
  Up_h_sum_ToT_mu.Write();
  Up_h_sum_ToT_casc_mu.Write();
  //------------------------------------------------------------
  Up_h_Ratio430.Write();
  Up_h_Ratio5.Write();
  Up_h_DiffDistance_mu.Write();
  Up_h_RatioCloseHits_mu.Write();
  Up_h_NtrackEarly.Write();
  Up_h_D.Write();
  Up_h_NtrackIT30.Write();
  Up_h_NtrackLate.Write();
  Up_h_myratio50_muon.Write();
  Up_h_ratio_cherenkov_lines.Write();
  Up_h_SnpeT.Write();
  Up_h_mean_tres_it.Write();
  Up_h_max_lik_up.Write();
  Up_h_max_lik_down.Write();
  Up_h_diff_theta.Write();
  Up_h_diff_dist_casc_mu.Write();
  Up_h_diff_dist_casc.Write();
  Up_h_ratio_closehits_cascmuon.Write();
  Up_h_ratio_closehits_casc.Write();
  Up_h_redToT_cascmuon.Write();
  Up_h_redToT_casc.Write();
  Up_h_myratio50_cascmuon.Write();
  Up_h_myratio50_casc.Write();
  Up_h_myratio30_cascmuon.Write();
  Up_h_myratio30_muon.Write();
  Up_h_myratio30_casc.Write();
  Up_h_min_diff_sollik.Write();
  Up_h_beta0_shower_deg.Write();
  Up_h_lik_shower.Write();
  Up_h_best_trk_pos_shower_z.Write();
  Up_h_normdlik.Write();
  Up_h_upsol.Write();
  Up_h_Slen.Write();
  Up_h_max_diff_sollik.Write();
  Up_h_diffangle_track_shower.Write();
  Up_h_myratio50_cascmuon_over_mu.Write();
  Up_h_myratio30_cascmuon_over_mu.Write();
  Up_h_myratio50_casc_over_mu.Write();
  Up_h_myratio30_casc_over_mu.Write();
  Up_h_ratio_closehits_cascmuon_over_mu.Write();
  Up_h_ratio_closehits_casc_over_mu.Write();
  Up_h_redToT_cascmuon_over_mu.Write();
  Up_h_redToT_casc_over_mu.Write();
  Up_h_bjorken_y.Write();
  Up_h_diffangle.Write();
  Up_h_LogE_mu.Write();
  Up_h_LogEbundle.Write();
  //Up_//h_logE_mu_max.Write();
  Up_h_cos_zen_mu.Write();
  //-------------------------------------------------------------

	
  Up_Lik_vs_NNhits.Write();//Nhits,jlik
  Up_Lik_vs_Q1value.Write();//Q1,jlik
  Up_NNhits_vs_Q1value.Write();//Q1,Nhits
  Up_Nlines_vs_zenith.Write();//zenith,nlines
  Up_Chere_Nlines_vs_zenith.Write();//zenith,chere_nlines
  Up_TrLengthIT_2_vs_Slen.Write();
  Up_Diffangle_vs_run.Write(); //run vs diffangle (at which run we have each diffangle)
  Up_LogEreco_vs_CosZen.Write(); 
  Up_LogEreco_vs_LogEreco2.Write();
  Up_LogEreco_vs_LogE_mu.Write();
  Up_LogEreco2_vs_LogE_mu.Write();
  Up_LogEreco_vs_LogE_neu.Write();
  Up_LogEreco2_vs_LogE_neu.Write();
  Up_LogEresolution_mu_vs_CosZen.Write();
  Up_LogEresolution_mu_cor_vs_CosZen.Write();
  Up_LogEresolution_neu_vs_CosZen.Write();
  Up_LogEresolution_neu_cor_vs_CosZen.Write();
     
      
  C_CosZen.Write();//cos_zen
  C_Lik.Write();//jlik
  C_LogBeta0.Write();//logbeta0
  C_Beta0_deg.Write();//logbeta0_deg
  C_DOMs.Write();//num_triggered_doms
  C_Cherenkov_DOMs.Write();//num_cherenkov_doms
  C_Q1value.Write();//Q1
  C_NNhits.Write();//Nhits
  C_LogEreco.Write();//log10(Ereco)
  C_TrLen.Write();//TrLengthIT_2
  C_Npe.Write();//SnpeT
  C_Ratiopmts.Write();//ratio_cherenkov_pmts
  C_Ratiodoms.Write();//ratio_cherenkov_doms
  C_Rvalue.Write();//R
  C_All_triggered_hits.Write();//num_triggered_hits
  C_Phi.Write();//fi
  C_Zenith_shower.Write();//zenith_shower
  C_DLik.Write();//dlik
  C_Ratio6.Write();//ratio6
  C_MaxZenSol.Write();//max_zen_sol
  C_RedToT_muon.Write();//redToT_muon
  C_Zenith.Write();//zenith
  C_Delta_ZeniTh_track_shower.Write();//delta_zenith_track_shower
  C_MinZenSol.Write();//min_zen_sol
  C_ItoLen.Write();//itoverlen
  C_Ratio330.Write();//ratio330
  C_Num_good_sol.Write();//num_of_good_sol
  C_PosZ.Write();//best_trk_pos_z
  C_LogEnu.Write();//logE
  C_LogAEnergy.Write();//logeReco_shower
  C_NLines.Write();//nlines
  C_LogEreco_cor.Write();
  C_Chere_NLines.Write();//chere_nlines
  C_Ratio_upper_hemisphere.Write();
  C_Ratio_lower_hemisphere.Write();
  C_Ratio_border_upper_hemisphere.Write();
  C_Ratio_border_lower_hemisphere.Write();
  C_Ratio1.Write();
  C_Ratio110.Write();
  C_Ratio130.Write();
  C_NtrackIT_ratio.Write();
  C_TrLenIT_3.Write();

  C_h_ToT_border_mu.Write();
  C_h_ToT_border_casc.Write();
  C_h_ToT_border_cascmu.Write();
  C_h_ToT_trig.Write();
  C_h_max_ToT_trig.Write();
  C_h_ToT_IT.Write();
  C_h_ToT_allIT.Write();
  C_h_sum_ToT_casc.Write();
  C_h_sum_ToT_mu.Write();
  C_h_sum_ToT_casc_mu.Write();
  //---------------------------------------------------------------
  C_h_Ratio430.Write();
  C_h_Ratio5.Write();
  C_h_DiffDistance_mu.Write();
  C_h_RatioCloseHits_mu.Write();
  C_h_NtrackEarly.Write();
  C_h_D.Write();
  C_h_NtrackIT30.Write();
  C_h_NtrackLate.Write();
  C_h_myratio50_muon.Write();
  C_h_ratio_cherenkov_lines.Write();
  C_h_SnpeT.Write();
  C_h_mean_tres_it.Write();
  C_h_max_lik_up.Write();
  C_h_max_lik_down.Write();
  C_h_diff_theta.Write();
  C_h_diff_dist_casc_mu.Write();
  C_h_diff_dist_casc.Write();
  C_h_ratio_closehits_cascmuon.Write();
  C_h_ratio_closehits_casc.Write();
  C_h_redToT_cascmuon.Write();
  C_h_redToT_casc.Write();
  C_h_myratio50_cascmuon.Write();
  C_h_myratio50_casc.Write();
  C_h_myratio30_cascmuon.Write();
  C_h_myratio30_muon.Write();
  C_h_myratio30_casc.Write();
  C_h_min_diff_sollik.Write();
  C_h_beta0_shower_deg.Write();
  C_h_lik_shower.Write();
  C_h_best_trk_pos_shower_z.Write();
  C_h_normdlik.Write();
  C_h_upsol.Write();
  C_h_Slen.Write();
  C_h_max_diff_sollik.Write();
  C_h_diffangle_track_shower.Write();
  C_h_myratio50_cascmuon_over_mu.Write();
  C_h_myratio30_cascmuon_over_mu.Write();
  C_h_myratio50_casc_over_mu.Write();
  C_h_myratio30_casc_over_mu.Write();
  C_h_ratio_closehits_cascmuon_over_mu.Write();
  C_h_ratio_closehits_casc_over_mu.Write();
  C_h_redToT_cascmuon_over_mu.Write();
  C_h_redToT_casc_over_mu.Write();
  C_h_bjorken_y.Write();
  C_h_diffangle.Write();
  C_h_LogE_mu.Write();
  C_h_LogEbundle.Write();
  //c_//h_logE_mu_max.Write();
  C_h_cos_zen_mu.Write();
  //-------------------------------------------------------------

  C_Lik_vs_NNhits.Write();//Nhits,jlik
  C_Lik_vs_Q1value.Write();//Q1,jlik
  C_NNhits_vs_Q1value.Write();//Q1,Nhits
  C_Nlines_vs_zenith.Write();//zenith,nlines
  C_Chere_Nlines_vs_zenith.Write();//zenith,chere_nlines
  C_TrLengthIT_2_vs_Slen.Write();
  C_Diffangle_vs_run.Write(); //run vs diffangle (at which run we have each diffangle)
  C_LogEreco_vs_CosZen.Write();
  C_LogEreco_vs_LogEreco2.Write();
  C_LogEreco_vs_LogE_mu.Write();
  C_LogEreco2_vs_LogE_mu.Write();
  C_LogEreco_vs_LogE_neu.Write();
  C_LogEreco2_vs_LogE_neu.Write();
  C_LogEresolution_mu_vs_CosZen.Write();
  C_LogEresolution_mu_cor_vs_CosZen.Write();
  C_LogEresolution_neu_vs_CosZen.Write();
  C_LogEresolution_neu_cor_vs_CosZen.Write();

      
  C_Up_CosZen.Write();//cos_zen
  C_Up_Lik.Write();//jlik
  C_Up_LogBeta0.Write();//logbeta0
  C_Up_Beta0_deg.Write();//logbeta0_deg
  C_Up_DOMs.Write();//num_triggered_doms
  C_Up_Cherenkov_DOMs.Write();//num_cherenkov_doms
  C_Up_Q1value.Write();//Q1
  C_Up_NNhits.Write();//Nhits
  C_Up_LogEreco.Write();//log10(Ereco)
  C_Up_TrLen.Write();//TrLengthIT_2
  C_Up_Npe.Write();//SnpeT
  C_Up_Ratiopmts.Write();//ratio_cherenkov_pmts
  C_Up_Ratiodoms.Write();//ratio_cherenkov_doms
  C_Up_Rvalue.Write();//R
  C_Up_All_triggered_hits.Write();//num_triggered_hits
  C_Up_Phi.Write();//fi
  C_Up_Zenith_shower.Write();//zenith_shower
  C_Up_DLik.Write();//dlik
  C_Up_Ratio6.Write();//ratio6
  C_Up_MaxZenSol.Write();//max_zen_sol
  C_Up_RedToT_muon.Write();//redToT_muon
  C_Up_Zenith.Write();//zenith
  C_Up_Delta_ZeniTh_track_shower.Write();//delta_zenith_track_shower
  C_Up_MinZenSol.Write();//min_zen_sol
  C_Up_ItoLen.Write();//itoverlen
  C_Up_Ratio330.Write();//ratio330
  C_Up_Num_good_sol.Write();//num_of_good_sol
  C_Up_PosZ.Write();//best_trk_pos_z
  C_Up_LogEnu.Write();//logE
  C_Up_LogAEnergy.Write();//logeReco_shower
  C_Up_NLines.Write();//nlines
  C_Up_LogEreco_cor.Write();
  C_Up_Chere_NLines.Write();//chere_nlines
  C_Up_Ratio_upper_hemisphere.Write();
  C_Up_Ratio_lower_hemisphere.Write();
  C_Up_Ratio_border_upper_hemisphere.Write();
  C_Up_Ratio_border_lower_hemisphere.Write();
  C_Up_Ratio1.Write();
  C_Up_Ratio110.Write();
  C_Up_Ratio130.Write();
  C_Up_NtrackIT_ratio.Write();
  C_Up_TrLenIT_3.Write();

  C_Up_h_ToT_border_mu.Write();
  C_Up_h_ToT_border_casc.Write();
  C_Up_h_ToT_border_cascmu.Write();
  C_Up_h_ToT_trig.Write();
  C_Up_h_max_ToT_trig.Write();
  C_Up_h_ToT_IT.Write();
  C_Up_h_ToT_allIT.Write();
  C_Up_h_sum_ToT_casc.Write();
  C_Up_h_sum_ToT_mu.Write();
  C_Up_h_sum_ToT_casc_mu.Write();
  //---------------------------------------------------------------
  C_Up_h_Ratio430.Write();
  C_Up_h_Ratio5.Write();
  C_Up_h_DiffDistance_mu.Write();
  C_Up_h_RatioCloseHits_mu.Write();
  C_Up_h_NtrackEarly.Write();
  C_Up_h_D.Write();
  C_Up_h_NtrackIT30.Write();
  C_Up_h_NtrackLate.Write();
  C_Up_h_myratio50_muon.Write();
  C_Up_h_ratio_cherenkov_lines.Write();
  C_Up_h_SnpeT.Write();
  C_Up_h_mean_tres_it.Write();
  C_Up_h_max_lik_up.Write();
  C_Up_h_max_lik_down.Write();
  C_Up_h_diff_theta.Write();
  C_Up_h_diff_dist_casc_mu.Write();
  C_Up_h_diff_dist_casc.Write();
  C_Up_h_ratio_closehits_cascmuon.Write();
  C_Up_h_ratio_closehits_casc.Write();
  C_Up_h_redToT_cascmuon.Write();
  C_Up_h_redToT_casc.Write();
  C_Up_h_myratio50_cascmuon.Write();
  C_Up_h_myratio50_casc.Write();
  C_Up_h_myratio30_cascmuon.Write();
  C_Up_h_myratio30_muon.Write();
  C_Up_h_myratio30_casc.Write();
  C_Up_h_min_diff_sollik.Write();
  C_Up_h_beta0_shower_deg.Write();
  C_Up_h_lik_shower.Write();
  C_Up_h_best_trk_pos_shower_z.Write();
  C_Up_h_normdlik.Write();
  C_Up_h_upsol.Write();
  C_Up_h_Slen.Write();
  C_Up_h_max_diff_sollik.Write();
  C_Up_h_diffangle_track_shower.Write();
  C_Up_h_myratio50_cascmuon_over_mu.Write();
  C_Up_h_myratio30_cascmuon_over_mu.Write();
  C_Up_h_myratio50_casc_over_mu.Write();
  C_Up_h_myratio30_casc_over_mu.Write();
  C_Up_h_ratio_closehits_cascmuon_over_mu.Write();
  C_Up_h_ratio_closehits_casc_over_mu.Write();
  C_Up_h_redToT_cascmuon_over_mu.Write();
  C_Up_h_redToT_casc_over_mu.Write();
  C_Up_h_bjorken_y.Write();
  C_Up_h_diffangle.Write();
  C_Up_h_LogE_mu.Write();
  C_Up_h_LogEbundle.Write();
  //C_Up_//h_logE_mu_max.Write();
  C_Up_h_cos_zen_mu.Write();
  //-------------------------------------------------------------

  C_Up_Lik_vs_NNhits.Write();//Nhits,jlik
  C_Up_Lik_vs_Q1value.Write();//Q1,jlik
  C_Up_NNhits_vs_Q1value.Write();//Q1,Nhits
  C_Up_Nlines_vs_zenith.Write();//zenith,nlines
  C_Up_Chere_Nlines_vs_zenith.Write();//zenith,chere_nlines
  C_Up_TrLengthIT_2_vs_Slen.Write();
  C_Up_Diffangle_vs_run.Write(); //run vs diffangle (at which run we have each diffangle)
  C_Up_LogEreco_vs_CosZen.Write();
  C_Up_LogEreco_vs_LogEreco2.Write();
  C_Up_LogEreco_vs_LogE_mu.Write();
  C_Up_LogEreco2_vs_LogE_mu.Write();
  C_Up_LogEreco_vs_LogE_neu.Write();
  C_Up_LogEreco2_vs_LogE_neu.Write();
  C_Up_LogEresolution_mu_vs_CosZen.Write();
  C_Up_LogEresolution_mu_cor_vs_CosZen.Write();
  C_Up_LogEresolution_neu_vs_CosZen.Write();
  C_Up_LogEresolution_neu_cor_vs_CosZen.Write();
      

  First_CosZen.Write();//cos_zen
  First_Lik.Write();//jlik
  First_LogBeta0.Write();//logbeta0
  First_Beta0_deg.Write();//logbeta0_deg
  First_DOMs.Write();//num_triggered_doms
  First_Cherenkov_DOMs.Write();//num_cherenkov_doms
  First_Q1value.Write();//Q1
  First_NNhits.Write();//Nhits
  First_LogEreco.Write();//log10(Ereco)    
  First_TrLen.Write();//TrLengthIT_2
  First_Npe.Write();//SnpeT
  First_Ratiopmts.Write();//ratio_cherenkov_pmts
  First_Ratiodoms.Write();//ratio_cherenkov_doms
  First_Rvalue.Write();//R
  First_All_triggered_hits.Write();//num_triggered_hits
  First_Phi.Write();//fi
  First_Zenith_shower.Write();//zenith_shower
  First_DLik.Write();//dlik
  First_Ratio6.Write();//ratio6
  First_MaxZenSol.Write();//max_zen_sol
  First_RedToT_muon.Write();//redToT_muon
  First_Zenith.Write();//zenith
  First_Delta_ZeniTh_track_shower.Write();//delta_zenith_track_shower
  First_MinZenSol.Write();//min_zen_sol
  First_ItoLen.Write();//itoverlen
  First_Ratio330.Write();//ratio330
  First_Num_good_sol.Write();//num_of_good_sol
  First_PosZ.Write();//best_trk_pos_z
  First_LogEnu.Write();//logE
  First_LogAEnergy.Write();//logeReco_shower
  First_NLines.Write();//nlines
  First_LogEreco_cor.Write();
  First_Chere_NLines.Write();//chere_nlines
  First_Ratio_upper_hemisphere.Write();
  First_Ratio_lower_hemisphere.Write();
  First_Ratio_border_upper_hemisphere.Write();
  First_Ratio_border_lower_hemisphere.Write();
  First_Ratio1.Write();
  First_Ratio110.Write();
  First_Ratio130.Write();
  First_NtrackIT_ratio.Write();
  First_TrLenIT_3.Write();

  First_h_ToT_border_mu.Write();
  First_h_ToT_border_casc.Write();
  First_h_ToT_border_cascmu.Write();
  First_h_ToT_trig.Write();
  First_h_max_ToT_trig.Write();
  First_h_ToT_IT.Write();
  First_h_ToT_allIT.Write();
  First_h_sum_ToT_casc.Write();
  First_h_sum_ToT_mu.Write();
  First_h_sum_ToT_casc_mu.Write();
  //---------------------------------------------------------------
  First_h_Ratio430.Write();
  First_h_Ratio5.Write();
  First_h_DiffDistance_mu.Write();
  First_h_RatioCloseHits_mu.Write();
  First_h_NtrackEarly.Write();
  First_h_D.Write();
  First_h_NtrackIT30.Write();
  First_h_NtrackLate.Write();
  First_h_myratio50_muon.Write();
  First_h_ratio_cherenkov_lines.Write();
  First_h_SnpeT.Write();
  First_h_mean_tres_it.Write();
  First_h_max_lik_up.Write();
  First_h_max_lik_down.Write();
  First_h_diff_theta.Write();
  First_h_diff_dist_casc_mu.Write();
  First_h_diff_dist_casc.Write();
  First_h_ratio_closehits_cascmuon.Write();
  First_h_ratio_closehits_casc.Write();
  First_h_redToT_cascmuon.Write();
  First_h_redToT_casc.Write();
  First_h_myratio50_cascmuon.Write();
  First_h_myratio50_casc.Write();
  First_h_myratio30_cascmuon.Write();
  First_h_myratio30_muon.Write();
  First_h_myratio30_casc.Write();
  First_h_min_diff_sollik.Write();
  First_h_beta0_shower_deg.Write();
  First_h_lik_shower.Write();
  First_h_best_trk_pos_shower_z.Write();
  First_h_normdlik.Write();
  First_h_upsol.Write();
  First_h_Slen.Write();
  First_h_max_diff_sollik.Write();
  First_h_diffangle_track_shower.Write();
  First_h_myratio50_cascmuon_over_mu.Write();
  First_h_myratio30_cascmuon_over_mu.Write();
  First_h_myratio50_casc_over_mu.Write();
  First_h_myratio30_casc_over_mu.Write();
  First_h_ratio_closehits_cascmuon_over_mu.Write();
  First_h_ratio_closehits_casc_over_mu.Write();
  First_h_redToT_cascmuon_over_mu.Write();
  First_h_redToT_casc_over_mu.Write();
  First_h_bjorken_y.Write();
  First_h_diffangle.Write();
  First_h_LogE_mu.Write();
  First_h_LogEbundle.Write();
  //First_//h_logE_mu_max.Write();
  First_h_cos_zen_mu.Write();
  //-------------------------------------------------------------
      
  First_Lik_vs_NNhits.Write();//Nhits,jlik
  First_Lik_vs_Q1value.Write();//Q1,jlik
  First_NNhits_vs_Q1value.Write();//Q1,Nhits
  First_Nlines_vs_zenith.Write();//zenith,nlines
  First_Chere_Nlines_vs_zenith.Write();//zenith,chere_nlines
  First_TrLengthIT_2_vs_Slen.Write();
  First_Diffangle_vs_run.Write(); //run vs diffangle (at which run we have each diffangle)
  First_LogEreco_vs_CosZen.Write();
  First_LogEreco_vs_LogEreco2.Write();
  First_LogEreco_vs_LogE_mu.Write();
  First_LogEreco2_vs_LogE_mu.Write();
  First_LogEreco_vs_LogE_neu.Write();
  First_LogEreco2_vs_LogE_neu.Write();
  First_LogEresolution_mu_vs_CosZen.Write();
  First_LogEresolution_mu_cor_vs_CosZen.Write();
  First_LogEresolution_neu_vs_CosZen.Write();
  First_LogEresolution_neu_cor_vs_CosZen.Write();

  First_Up_CosZen.Write();//cos_zen
  First_Up_Lik.Write();//jlik
  First_Up_LogBeta0.Write();//logbeta0
  First_Up_Beta0_deg.Write();//logbeta0_deg
  First_Up_DOMs.Write();//num_triggered_doms
  First_Up_Cherenkov_DOMs.Write();//num_cherenkov_doms
  First_Up_Q1value.Write();//Q1
  First_Up_NNhits.Write();//Nhits
  First_Up_LogEreco.Write();//log10(Ereco)
  First_Up_TrLen.Write();//TrLengthIT_2
  First_Up_Npe.Write();//SnpeT
  First_Up_Ratiopmts.Write();//ratio_cherenkov_pmts
  First_Up_Ratiodoms.Write();//ratio_cherenkov_doms
  First_Up_Rvalue.Write();//R
  First_Up_All_triggered_hits.Write();//num_triggered_hits
  First_Up_Phi.Write();//fi
  First_Up_Zenith_shower.Write();//zenith_shower
  First_Up_DLik.Write();//dlik
  First_Up_Ratio6.Write();//ratio6
  First_Up_MaxZenSol.Write();//max_zen_sol
  First_Up_RedToT_muon.Write();//redToT_muon
  First_Up_Zenith.Write();//zenith
  First_Up_Delta_ZeniTh_track_shower.Write();//delta_zenith_track_shower
  First_Up_MinZenSol.Write();//min_zen_sol
  First_Up_ItoLen.Write();//itoverlen
  First_Up_Ratio330.Write();//ratio330
  First_Up_Num_good_sol.Write();//num_of_good_sol
  First_Up_PosZ.Write();//best_trk_pos_z
  First_Up_LogEnu.Write();//logE
  First_Up_LogAEnergy.Write();//logeReco_shower
  First_Up_NLines.Write();//nlines
  First_Up_LogEreco_cor.Write();
  First_Up_Chere_NLines.Write();//chere_nlines
  First_Up_Ratio_upper_hemisphere.Write();
  First_Up_Ratio_lower_hemisphere.Write();
  First_Up_Ratio_border_upper_hemisphere.Write();
  First_Up_Ratio_border_lower_hemisphere.Write();
  First_Up_Ratio1.Write();
  First_Up_Ratio110.Write();
  First_Up_Ratio130.Write();
  First_Up_NtrackIT_ratio.Write();
  First_Up_TrLenIT_3.Write();

  First_Up_h_ToT_border_mu.Write();
  First_Up_h_ToT_border_casc.Write();
  First_Up_h_ToT_border_cascmu.Write();
  First_Up_h_ToT_trig.Write();
  First_Up_h_max_ToT_trig.Write();
  First_Up_h_ToT_IT.Write();
  First_Up_h_ToT_allIT.Write();
  First_Up_h_sum_ToT_casc.Write();
  First_Up_h_sum_ToT_mu.Write();
  First_Up_h_sum_ToT_casc_mu.Write();
  //---------------------------------------------------------------
  First_Up_h_Ratio430.Write();
  First_Up_h_Ratio5.Write();
  First_Up_h_DiffDistance_mu.Write();
  First_Up_h_RatioCloseHits_mu.Write();
  First_Up_h_NtrackEarly.Write();
  First_Up_h_D.Write();
  First_Up_h_NtrackIT30.Write();
  First_Up_h_NtrackLate.Write();
  First_Up_h_myratio50_muon.Write();
  First_Up_h_ratio_cherenkov_lines.Write();
  First_Up_h_SnpeT.Write();
  First_Up_h_mean_tres_it.Write();
  First_Up_h_max_lik_up.Write();
  First_Up_h_max_lik_down.Write();
  First_Up_h_diff_theta.Write();
  First_Up_h_diff_dist_casc_mu.Write();
  First_Up_h_diff_dist_casc.Write();
  First_Up_h_ratio_closehits_cascmuon.Write();
  First_Up_h_ratio_closehits_casc.Write();
  First_Up_h_redToT_cascmuon.Write();
  First_Up_h_redToT_casc.Write();
  First_Up_h_myratio50_cascmuon.Write();
  First_Up_h_myratio50_casc.Write();
  First_Up_h_myratio30_cascmuon.Write();
  First_Up_h_myratio30_muon.Write();
  First_Up_h_myratio30_casc.Write();
  First_Up_h_min_diff_sollik.Write();
  First_Up_h_beta0_shower_deg.Write();
  First_Up_h_lik_shower.Write();
  First_Up_h_best_trk_pos_shower_z.Write();
  First_Up_h_normdlik.Write();
  First_Up_h_upsol.Write();
  First_Up_h_Slen.Write();
  First_Up_h_max_diff_sollik.Write();
  First_Up_h_diffangle_track_shower.Write();
  First_Up_h_myratio50_cascmuon_over_mu.Write();
  First_Up_h_myratio30_cascmuon_over_mu.Write();
  First_Up_h_myratio50_casc_over_mu.Write();
  First_Up_h_myratio30_casc_over_mu.Write();
  First_Up_h_ratio_closehits_cascmuon_over_mu.Write();
  First_Up_h_ratio_closehits_casc_over_mu.Write();
  First_Up_h_redToT_cascmuon_over_mu.Write();
  First_Up_h_redToT_casc_over_mu.Write();
  First_Up_h_bjorken_y.Write();
  First_Up_h_diffangle.Write();
  First_Up_h_LogE_mu.Write();
  First_Up_h_LogEbundle.Write();
  //First_Up_//h_logE_mu_max.Write();
  First_Up_h_cos_zen_mu.Write();
  //-------------------------------------------------------------

  First_Up_Lik_vs_NNhits.Write();//Nhits,jlik
  First_Up_Lik_vs_Q1value.Write();//Q1,jlik
  First_Up_NNhits_vs_Q1value.Write();//Q1,Nhits
  First_Up_Nlines_vs_zenith.Write();//zenith,nlines
  First_Up_Chere_Nlines_vs_zenith.Write();//zenith,chere_nlines
  First_Up_TrLengthIT_2_vs_Slen.Write();
  First_Up_Diffangle_vs_run.Write(); //run vs diffangle (at which run we have each diffangle)
  First_Up_LogEreco_vs_CosZen.Write();
  First_Up_LogEreco_vs_LogEreco2.Write();
  First_Up_LogEreco_vs_LogE_mu.Write();
  First_Up_LogEreco2_vs_LogE_mu.Write();
  First_Up_LogEreco_vs_LogE_neu.Write();
  First_Up_LogEreco2_vs_LogE_neu.Write();
  First_Up_LogEresolution_mu_vs_CosZen.Write();
  First_Up_LogEresolution_mu_cor_vs_CosZen.Write();
  First_Up_LogEresolution_neu_vs_CosZen.Write();
  First_Up_LogEresolution_neu_cor_vs_CosZen.Write();

  //comon level
  Common_CosZen.Write();
  Common_Lik.Write();
  Common_LogBeta0.Write();
  Common_LogEreco.Write();
  Common_TrLen.Write();
  Common_LogEresolution_mu_vs_CosZen.Write();
  Common_LogEresolution_neu_vs_CosZen.Write();
  //---------------- 

  Rate3Dmuon.Write();
  C_Rate3Dmuon.Write();
  First_Rate3Dmuon.Write();
  Up_Rate3Dmuon.Write();
  C_Up_Rate3Dmuon.Write();
  First_Up_Rate3Dmuon.Write();
  Rate3Dshower.Write();
  C_Rate3Dshower.Write();
  First_Rate3Dshower.Write();
  Up_Rate3Dshower.Write();
  C_Up_Rate3Dshower.Write();
  First_Up_Rate3Dshower.Write();
  RateMXshower.Write();
  C_RateMXshower.Write();
  First_RateMXshower.Write();
  Up_RateMXshower.Write();
  C_Up_RateMXshower.Write();
  First_Up_RateMXshower.Write();

  Allreco_rate.Write();
  C_Allreco_rate.Write();
  First_Allreco_rate.Write();

  Up_Allreco_rate.Write();
  C_Up_Allreco_rate.Write();
  First_Up_Allreco_rate.Write();

  
 
  //---------------------------------------

  

  std::cout<<"-- Everything's fine! --"<<std::endl;

  f_out->Close();




}//end of main

float icecube_flux_diff2019(float E) {
  float E_over_100_TeV = E / 100000.0;
  float f = 1.44 * std::pow(10, -18) * std::pow(E_over_100_TeV, -2.28);  // flux  from IceCube diffuse astrophysics numuCC analysis in GeV^-1 s^-1 cm^-2 sr^-1
  return 0.5 * 10000.0 * f;                                              // have flux in GeV^-1 s^-1 m^-2 sr^-1
}

float icecube_flux_diff2021(float E){
  float E_over_100_TeV = E / 100000.0;
  float f = 1.44 * std::pow(10, -18) * std::pow(E_over_100_TeV, -2.37);  // flux  from IceCube diffuse astrophysics numuCC analysis in GeV^-1 s^-1 cm^-2 sr^-1
  return 0.5 * 10000.0 * f;                                              // have flux in GeV^-1 s^-1 m^-2 sr^-1
}

