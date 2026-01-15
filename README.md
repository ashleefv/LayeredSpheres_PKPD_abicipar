# LayeredSpheres_PKPD_abicipar
## Overview
This project combines a mathematical model of drug release from single‑layer and bi‑layer polymeric microsphere DDSs with a pharmacokinetic/pharmacodynamic (PK/PD) model of abicipar pegol and VEGF inhibition. The model tracks abicipar distribution across ocular tissues (vitreous, aqueous humor, retina, choroid) and serum, and evaluates VEGF suppression in key compartments.
## Authors
Sarita Das<sup>a</sup>, Mohammad Aminul Islam<sup>f</sup>, Yaman Oklla<sup>a</sup>, Koki Kanehira<sup>a</sup>,
Eduardo A. Chacin Ruiz<sup>a</sup>, Katelyn E. Swindle-Reilly<sup>e,f,g</sup>, Ashlee N. Ford Versypt<sup>a,b,c,d</sup><br/>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, 14260, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, 14260, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, 14260, USA<br/>
<sup>d</sup>Witebsky Center for Microbial Pathogenesis and Immunology, University at Buffalo, The State University of New York, Buffalo, NY, 14203, USA
<sup>e</sup>Department of Pharmaceutical Sciences, University at Buffalo, The State University of New York, Buffalo, NY, 14215, USA<br/>
<sup>f</sup>Department of Biomedical Engineering, University of Delaware, Newark, DE, 19716, USA
<sup>g</sup>William G. Lowrie Department of Chemical and Biomolecular Engineering, The Ohio State University, Columbus, OH, 43210, USA<br/>
<sup>h</sup>Department of Ophthalmology and Visual Sciences, The Ohio State University, Columbus, OH, 43210, USA<br/>
<sup>i</sup>Department of Biomedical Engineering, The Ohio State University, Columbus, OH, 43210, USA<br/>

## Manuscript

## Scripts

* **without_DDS.m** This file simulates abicipar pharmacokinetics and computes the time required for VEGF levels to return to 10% and 50% of their baseline concentrations following various abicipar doses administered without a DDS.

* **DDS_doses.m** This file simulates abicipar pharmacokinetics and computes the time required for VEGF levels to return to 10% and 50% of their baseline concentrations following various abicipar doses administered with DDS.

* **chitosan_single.m** This file simulates abicipar pharmacokinetics and the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 0.1 mg from single-layered chitosan DDS. The radius of the chitosan layer is varied in the simulation.

* **PCL_single.m** This file simulates abicipar pharmacokinetics and the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 0.1 mg from single-layered PCL DDS. The radius of the PCL layer is varied in the simulation.

* **bi_layer_chitosan.m** This file simulates abicipar pharmacokinetics and the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 0.1 mg from bi-layered chitosan-PCL DDS. The radius of the chitosan layer is varied in the simulation.

* **bi_layer_PCL.m** This file simulates abicipar pharmacokinetics and the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 0.1 mg from bi-layered chitosan-PCL DDS. The thiskness of the PCL layer is varied in the simulation.

* **bi_layer_changing_both.m** This file simulates abicipar pharmacokinetics and the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 0.1 mg from bi-layered chitosan-PCL DDS. The radius of the chitosan layer and the thiskness of the PCL layer is varied simultaneously in the simulation.

* **FD_spheres_variable_diffusivity_two_spheres.m** This is the function file to solve the PDE of bi-layered chitosan-PCL core-shell DDS for Fickian diffusion within a radially symmetric sphere

* **solve_FD_spheres_variable_diffusivity.m** This file solves the PDE of bi-layered chitosan-PCL core-shell DDS for Fickian diffusion within a radially symmetric sphere

* **plot_all.m** This file generates manuscrpt figures 2-9.

## Acknowledgements
This work was supported by National Institutes of Health grant R35GM133763 to ANFV, R01EB032870 to KESR and ANFV, Owen Locke Foundation to KESR, and the University at Buffalo.