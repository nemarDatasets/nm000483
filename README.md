[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000483-blue)](https://doi.org/10.82901/nemar.nm000483)

# sEEG seizure-onset and interictal segments of one polymicrogyria patient, with electrode coordinates (Vakilna et al., 2025) - derivative

Release: Zenodo 14807262 "Time-Frequency Fingerprint Analysis in SEEG Source-Space to Identify the Epileptogenic Zone" (Vakilna, Yash Shashank, Atilgan, Deniz, Hampson, Johnson P, Chinara, Chinmay, Medani, Takfarinas, Leahy, Richard, Lacuey Lecumberri, Nuria, Lhatoo, Samden, Pati, Sandipan, Mosher, John, Gavvala, Jay; doi:10.5281/zenodo.14807262, CC-BY-4.0), the Brainstorm database
(`BST_dataset.zip`) behind the case study doi:10.1002/acn3.70115. Zenodo description: "This case study highlights the application of seizure fingerprint analysis in the source-space of SEEG data to accurately localize the epileptogenic zone (EZ) in patients with complex cortical malformations. A 25-year-old female with extensive bilateral perisylvian polymicrogyria presented with intractable focal seizures. The source-level analysis performed in Brainstorm using the sLORETA imaging algorithm subsequently showed EZ fingerprint analysis. The patient underwent MR-guided laser interstitial thermal therapy, targeting the identified EZ, resulting in postoperative seizure freedom with minimal complications. Extending the analysis from the sensor-space to the source-space could further enhance surgical planning and improve outcomes in complex epilepsy cases."

## Content (1 participant, 5 recordings, 1.2 min of sensor-space data)
The release contains three short sensor-space segments of the clinical stereo-EEG (Nihon Kohden EEG-1200A, 2000 Hz) that
the authors imported into Brainstorm, notch-filtered (180-900 Hz harmonics) and resampled to 1000 Hz, and their bipolar
('bipolar 2') montage versions. They are written here as BrainVision float32 in microvolts (Brainstorm stores volts),
`acq-monopolar` (267 channels: 228 located sEEG contacts, 20 sEEG contacts without location, scalp EEG, ECG, DC and
marker channels, as typed in the Brainstorm channel file) and `acq-bipolar` (206 bipolar pairs). This is a DERIVATIVE
dataset; the full raw recording is not in the release. Brainstorm events (clinical marks, the authors' onset/wave/spike
markers and channel annotations) are in events.tsv; the time axis of each segment (0 = reference event) is kept in the
`brainstorm_time` column.
  - `sub-01_task-seizureonsetlvfa_acq-monopolar`: LVFA_and_wave/data_sEEG_onset__SPS_PLON_trial001_notch_resample.mat (sEEG onset: SPS/PLON (#1) | notch(180Hz 300Hz 420Hz 540Hz 660Hz 780Hz 900Hz)), 267 ch, 1000 Hz, 30 s, time axis -15..15 s, 7 events
  - `sub-01_task-seizureonsetlvfa_acq-bipolar`: LVFA_and_wave_bipolar_2/data_sEEG_onset__SPS_PLON_trial001_notch_resample_montage.mat (sEEG onset: SPS/PLON (#1) | notch(180Hz 300Hz 420Hz 540Hz 660Hz 780Hz 900Hz) | bipolar 2), 206 ch, 1000 Hz, 30 s, time axis -15..15 s, 7 events
  - `sub-01_task-seizureonsetspike_acq-monopolar`: ictal_repetitive_spike/data_EEG_onset_trial001_notch_resample.mat (sEEG onset (#1) | notch(180Hz 300Hz 420Hz 540Hz 660Hz 780Hz 900Hz)), 267 ch, 1000 Hz, 30 s, time axis -15..15 s, 6 events
  - `sub-01_task-seizureonsetspike_acq-bipolar`: ictal_repetitive_spike_bipolar_2/data_EEG_onset_trial001_notch_resample_montage.mat (sEEG onset (#1) | notch(180Hz 300Hz 420Hz 540Hz 660Hz 780Hz 900Hz) | bipolar 2), 206 ch, 1000 Hz, 30 s, time axis -15..15 s, 6 events
  - `sub-01_task-interictalspike_acq-monopolar`: interictal_spike/data_Spike_trial001_notch_resample.mat (Spike (#1) | notch(180Hz 300Hz 420Hz 540Hz 660Hz 780Hz 900Hz)), 267 ch, 1000 Hz, 10 s, time axis -5..5 s, 3 events

## Localisation
`electrodes.tsv` carries the Brainstorm channel locations (subject coordinate system SCS, native space, mm) of 228 sEEG
contacts of 16 depth electrodes (AH, ALAT, ALON, MC, MPM, MPS, PH, PIN, PLAT, PLON, POP, PSIN, SMA, SPL, SPM, SPS); electrode geometry (contact number,
spacing, diameter) is in `sourcedata/zenodo-14807262/intra_electrodes.json`. The subject's T1 MRI, head/cortex surfaces and
volume atlases of the release are NOT redistributed (the MRI is not defaced); they, and the authors' source-space results
(OpenMEEG head model, sLORETA kernels, time-frequency maps, scout time series), remain in the Zenodo record.

## Provenance
Converted by the ieeg-nemar-20261005 campaign (lane o5c-B) with `code/o5c_B_conv_lit187.py`; round-trip checks (MNE
read-back) in code/conversion_report.json. sourcedata holds the converted Brainstorm data/channel/study files verbatim.

## RecordingType (18 Oct 2026)

5 data files were marked `RecordingType` epoched without a usable epoch definition (no `EpochLength`, or `EpochLength` equal to the whole file). Each of them is ONE contiguous excerpt (no segment boundary in the .vmrk, no gap): there is nothing to cut into epochs, `read_raw_bids` refuses epoched files and `read_epochs_bids` failed (missing `EpochLength`, or an epoch dropped for every event that does not start at sample 0). `RecordingType` is now continuous (and `EpochLength` removed where it was present, 5 files); the files read with `read_raw_bids` and every event stays an annotation. The excerpts are still the authors' fixed windows (see above); only the BIDS label changed.

## Channel types and non-neural channel values (fix-20261030, Oct 2026)

- `sub-01_task-seizureonsetspike_acq-monopolar`: the 11 scalp electrodes `F3`, `C3`, `P3`, `O1`, `F4`, `C4`, `P4`, `O2`,
  `FZ`, `CZ`, `PZ` were typed MISC, because that file's Brainstorm channel file (`ictal_repetitive_spike/channel_02.mat`)
  types them MISC. The channel files of the other two monopolar segments of the same recording (`LVFA_and_wave/channel.mat`,
  `interictal_spike/channel.mat`) type the same 11 electrodes EEG, so they are now typed EEG here too.
- Values of the non-neural channels: all channels were converted from Brainstorm's volts to microvolts (x 1e6), including
  the channels that do not carry voltages. The marker channel `Events` (TRIG) therefore holds values up to 3.5e10, and the
  MISC channels (DC inputs, blood pressure, CO2, SpO2, respiration and other monitors) up to 4.1e5 (maximum |value| per
  monopolar file: TRIG 3.50e10; MISC 6.1e4, 4.1e5 and 1.4e5; SEEG at most 3,335; check of 30 Oct 2026). These channels
  have `units` n/a; their values are not microvolts and should not be used as signals.
