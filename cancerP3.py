#!/usr/bin/env python
# coding: utf-8

# In[48]:


from pyopenms import *
import matplotlib.pyplot as plt

protein_ids = []
peptide_ids = []
SimpleSearchEngineAlgorithm().search("cc3.mzML", "cancer.fasta", protein_ids, peptide_ids)
for peptide_id in peptide_ids:
    print ("___________________________________________________________________________")
    print ("Peptide Index m/z:", peptide_id.getMetaValue("scan_index"))
    for hit in peptide_id.getHits():
        print(" - Peptide hit sequence:", hit.getSequence())
        print(" - Peptide hit rank:", hit.getRank())
        string= hit.getSequence().toString()
        tsg = TheoreticalSpectrumGenerator()
        theo_spec = MSSpectrum()
        p = Param()
        p.setValue("add_y_ions", "true")
        p.setValue("add_b_ions", "true")
        p.setValue("add_metainfo", "true")
        tsg.setParameters(p)
        peptide = AASequence.fromString(hit.getSequence().toString())
        tsg.getSpectrum(theo_spec, peptide, 1, 2)
        # Iterate over annotated ions and their masses
        print("Spectrum 1 of", peptide, "has", theo_spec.size(), "peaks.")
        for ion, peak in zip(theo_spec.getStringDataArrays()[0], theo_spec):
            print(ion.decode(), "is generated at m/z", peak.getMZ())
        exp = MSExperiment()
        MzMLFile().load("cc3.mzML", exp)
        spectra = exp.getSpectrum(peptide_id.getMetaValue("scan_index"))
        alignment = []
        spa = SpectrumAlignment()
        p = spa.getParameters()
        p.setValue("tolerance", 0.5)
        p.setValue("is_relative_tolerance", "false")
        spa.setParameters(p)
        # align both spectra
        spa.getSpectrumAlignment(alignment, theo_spec, spectra)

        # Print matching ions and mz from theoretical spectrum
        print("Number of matched peaks: " + str(len(alignment)))
        print("ion\ttheo. m/z\tobserved m/z")

        for theo_idx, obs_idx in alignment:
            ion_name = theo_spec.getStringDataArrays()[0][theo_idx].decode()
            ion_charge = theo_spec.getIntegerDataArrays()[0][theo_idx]
            print(ion_name + "\t" + str(ion_charge) + "\t"
                  + str(theo_spec[theo_idx].getMZ())
                  + "\t" + str(spectra[obs_idx].getMZ()))

        theo_mz, theo_int, obs_mz, obs_int = [], [], [], []
        for theo_idx, obs_idx in alignment:
            theo_mz.append(theo_spec[theo_idx].getMZ())
            theo_int.append(theo_spec[theo_idx].getIntensity())
            obs_mz.append(spectra[obs_idx].getMZ())
            obs_int.append(spectra[obs_idx].getIntensity())
        title = f'{string},{peptide_id.getMetaValue("scan_index")}'
        flag=False
        for i in obs_int:
            if(i>0):
                flag=True
                break
        if(flag==True):
           mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title)
        else:
           flag=False

def mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title):
  obs_int = [element / max(obs_int) for element in obs_int]  # relative intenstiy
  theo_int = [element * -1 for element in theo_int]  # invert the intensity for the mirror plot
  plt.figure(figsize=(12, 8))
  plt.bar(obs_mz, obs_int, width=3.0)
  plt.bar(theo_mz, theo_int, width=3.0)
  plt.title(title)
  plt.ylabel('intensity')
  plt.xlabel('m/z')
  plt.show()        


# In[ ]:





# In[ ]:




