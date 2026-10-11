"""o5c lane B: src-literature-187 (Zenodo 14807262, Vakilna et al. 2025, 'Time-Frequency Fingerprint Analysis in SEEG
Source-Space to Identify the Epileptogenic Zone') -> iEEG-BIDS DERIVATIVE. The release is a Brainstorm database
(BST_dataset.zip) of one patient: three sEEG segments imported from Nihon Kohden/EDF, notch-filtered and resampled to
1000 Hz by the authors (data_*_notch_resample.mat, 267 channels x 30/30/10 s), their bipolar-montage versions, a channel
file with SCS coordinates (Brainstorm subject coordinate system, metres) for 228 contacts and 16 IntraElectrodes
definitions, plus source-space results (head model, sLORETA kernels, time-frequency maps, scouts) and the subject's
T1 MRI / surfaces. Converted: the three sensor-space segments (acq-monopolar) and their bipolar versions (acq-bipolar) as
BrainVision float32 in uV; Brainstorm Events -> events.tsv; Channel.Loc -> electrodes.tsv (space-Other = Brainstorm SCS,
mm). Not included: the MRI and surfaces (not defaced), source-space results (in the Zenodo record).
Usage: python o5c_B_conv_lit187.py src-literature-187"""
import glob, json, os, re, shutil, sys, urllib.request
import numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_B_common import *  # noqa
from b3w10_common import write_run  # noqa

CID = sys.argv[1] if len(sys.argv) > 1 else "src-literature-187"
SRC = src_dir(CID); OUT = fresh_out(CID); REP = rep_dir(CID)
B = os.path.join(SRC, "_x", "BST_dataset.zip.d")
ZEN = "14807262"; DOI = f"10.5281/zenodo.{ZEN}"; URL = f"https://zenodo.org/records/{ZEN}"
CONV = "harness/o5c_B_conv_lit187.py (ieeg-nemar-20261005 campaign, lane o5c-B)"
try:
    rec = json.load(urllib.request.urlopen(f"https://zenodo.org/api/records/{ZEN}", timeout=60))
except Exception as e:
    rec = {"error": repr(e)}
SD = os.path.join(OUT, "sourcedata", f"zenodo-{ZEN}"); os.makedirs(SD, exist_ok=True)
wjson(os.path.join(SD, f"zenodo_record_{ZEN}.json"), rec)
meta = rec.get("metadata", {})
authors = [c["name"] for c in meta.get("creators", [])]
rep = {"cand_id": CID, "converter": CONV, "licence": {"value": (meta.get("license") or {}).get("id", "cc-by-4.0"), "read_from": "zenodo api metadata.license.id"}, "runs": [], "sourcedata": []}
sub = "sub-01"; D = os.path.join(OUT, sub, "ieeg"); os.makedirs(D, exist_ok=True)
SEG = [("LVFA_and_wave", "data_sEEG_onset__SPS_PLON_trial001_notch_resample.mat", "LVFA_and_wave_bipolar_2", "data_sEEG_onset__SPS_PLON_trial001_notch_resample_montage.mat", "seizureonsetlvfa", "Seizure onset with low-voltage fast activity and 'wave' (release study folder LVFA_and_wave; segment -15..+15 s around the sEEG onset marked SPS/PLON)"),
       ("ictal_repetitive_spike", "data_EEG_onset_trial001_notch_resample.mat", "ictal_repetitive_spike_bipolar_2", "data_EEG_onset_trial001_notch_resample_montage.mat", "seizureonsetspike", "Seizure onset with repetitive spiking (release study folder ictal_repetitive_spike; -15..+15 s around the sEEG onset)"),
       ("interictal_spike", "data_Spike_trial001_notch_resample.mat", None, None, "interictalspike", "Interictal spike (release study folder interictal_spike; -5..+5 s around the spike)")]
TYPEMAP = {"SEEG": "SEEG", "SEEG_NO_LOC": "SEEG", "EEG": "EEG", "EEG REF": "REF", "ECG": "ECG", "DC": "MISC", "MISC": "MISC", "STIM": "TRIG"}


def load_channel(p):
    m = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
    out = []
    for c in np.atleast_1d(m["Channel"]):
        loc = np.asarray(c.Loc, dtype=float).flatten()
        out.append({"name": str(c.Name), "bst_type": str(c.Type), "group": str(c.Group) if np.asarray(c.Group).size else "", "loc_m": loc[:3].tolist() if loc.size >= 3 else None})
    ie = []
    for x in np.atleast_1d(m.get("IntraElectrodes", [])):
        if hasattr(x, "_fieldnames"):
            ie.append({"Name": str(x.Name), "Type": str(x.Type), "ContactNumber": float(np.asarray(x.ContactNumber).flatten()[0]) if np.asarray(x.ContactNumber).size else None,
                       "ContactSpacing_m": float(np.asarray(x.ContactSpacing).flatten()[0]) if np.asarray(x.ContactSpacing).size else None,
                       "ContactDiameter_m": float(np.asarray(x.ContactDiameter).flatten()[0]) if np.asarray(x.ContactDiameter).size else None,
                       "ContactLength_m": float(np.asarray(x.ContactLength).flatten()[0]) if np.asarray(x.ContactLength).size else None,
                       "ElecLength_m": float(np.asarray(x.ElecLength).flatten()[0]) if np.asarray(x.ElecLength).size else None})
    return out, ie


elec_rows = {}; ie_all = None
for folder, fn, bfolder, bfn, task, tdesc in SEG:
    for acq, fo, f_ in (("monopolar", folder, fn), ("bipolar", bfolder, bfn)):
        if fo is None:
            continue
        p = os.path.join(B, "data", "Subject01", fo, f_)
        chp = sorted(glob.glob(os.path.join(B, "data", "Subject01", fo, "channel*.mat")))[0]
        chans, ie = load_channel(chp)
        if ie_all is None and acq == "monopolar":
            ie_all = ie
        m = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
        F = np.asarray(m["F"], dtype=np.float64); t = np.asarray(m["Time"]).flatten()
        assert F.shape[0] == len(chans) == len(np.asarray(m["ChannelFlag"]).flatten()), (F.shape, len(chans))
        fs = float(round(1.0 / float(np.median(np.diff(t))), 6))
        flag = np.asarray(m["ChannelFlag"]).flatten()
        names = [c["name"] for c in chans]
        types = [TYPEMAP.get(c["bst_type"], "MISC") for c in chans]
        units = ["µV" if ty in ("SEEG", "EEG", "ECG", "REF") else "n/a" for ty in types]
        X = (F * 1e6).T  # Brainstorm stores volts
        ch_extra = [{"status": "bad" if flag[i] < 0 else "good", "status_description": "ChannelFlag -1 in the Brainstorm file" if flag[i] < 0 else "n/a",
                     "group": chans[i]["group"] or "n/a", "reference": ("bipolar (adjacent contacts, Brainstorm 'bipolar 2' montage)" if acq == "bipolar" else "n/a"),
                     "description": f"Brainstorm channel type {chans[i]['bst_type']}" + ("; no location in the channel file" if chans[i]["bst_type"] == "SEEG_NO_LOC" else "")} for i in range(len(chans))]
        ev = []
        for e in np.atleast_1d(m.get("Events", [])):
            if not hasattr(e, "label"):
                continue
            tm = np.atleast_2d(np.asarray(e.times, dtype=float))
            if tm.shape[0] == 2:  # extended events: onset/offset rows
                for a, b_ in zip(tm[0], tm[1]):
                    ev.append([round(float(a - t[0]), 6), round(float(b_ - a), 6), int(round((a - t[0]) * fs)), str(e.label), round(float(a), 6)])
            else:
                for a in tm.flatten():
                    ev.append([round(float(a - t[0]), 6), "n/a", int(round((a - t[0]) * fs)), str(e.label), round(float(a), 6)])
        ev.sort(key=lambda r: r[0])
        stem = f"{sub}_task-{task}_acq-{acq}"
        comment = str(m.get("Comment", ""))
        hist = [[str(x) for x in np.atleast_1d(row)] for row in np.atleast_2d(m.get("History", []))]
        proc = [h[2] for h in hist if len(h) >= 3 and h[1] in ("process", "import")]
        r = write_run(D, stem, X, fs, names, units, types, {
            "TaskName": task, "TaskDescription": tdesc + ". No task: clinical sEEG monitoring segment selected by the authors.",
            "PowerLineFrequency": 60, "SoftwareFilters": {"Brainstorm (authors)": {"steps": [re.sub(r"[A-Z]:\\\\\S+", "<path>", s_) for s_ in proc][:12], "comment": re.sub(r"[A-Z]:\\\\\S+", "<path>", comment)}},
            "HardwareFilters": "n/a", "iEEGReference": "bipolar (Brainstorm 'bipolar 2' montage)" if acq == "bipolar" else "n/a (recording reference of the clinical system; channel typed EEG REF is the reference electrode)",
            "iEEGPlacementScheme": f"Stereo-EEG, 16 depth electrodes ({', '.join(x['Name'] for x in (ie_all or []))}); release case: 25-year-old female with bilateral perisylvian polymicrogyria",
            "Manufacturer": "Nihon Kohden", "ManufacturersModelName": str(m.get("Device", "n/a")), "ElectricalStimulation": False, "InstitutionName": "n/a",
            "RecordingType": "epoched", "EpochLength": round(float(t[-1] - t[0] + 1 / fs), 3)}, ch_extra=ch_extra,
            comment=f"{CONV}: {fo}/{f_}", markers=[("Stimulus", f"t0 {t[0]:.3f} s relative to the Brainstorm event", int(round(-t[0] * fs)), 1)] if t[0] < 0 else None)
        wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "sample", "trial_type", "brainstorm_time"], ev)
        wjson(os.path.join(D, stem + "_events.json"), {"trial_type": {"Description": "Brainstorm event label as released (MARK ON/OFF = clinical system marks; 'wave', 'spiking', 'spike' = authors' markers; 'PIN 1-6', 'PLON 1/2', 'PLAT 6-8' = authors' channel annotations at onset; transient_notch = filter edge artefact window)"},
              "brainstorm_time": {"Description": "event time on the Brainstorm time axis (0 = the segment's reference event, e.g. seizure onset)", "Units": "s"},
              "duration": {"Description": "n/a for point events; extended events carry their length"}})
        r.update(sub=sub, task=task, acq=acq, stem=stem, source=f"{fo}/{f_}", n_events=len(ev), time_axis=[float(t[0]), float(t[-1])], n_bad=int((flag < 0).sum()), comment=comment)
        rep["runs"].append(r)
        print(stem, json.dumps({k: r[k] for k in ("n_samples", "sfreq", "n_channels", "roundtrip_ok", "n_events")}, default=str), flush=True)
        if acq == "monopolar":
            for c, ty in zip(chans, types):
                if ty in ("SEEG", "EEG", "ECG", "REF") and c["name"] not in elec_rows:
                    loc = c["loc_m"]
                    has = loc is not None and len(loc) == 3 and any(abs(v) > 0 for v in loc)
                    elec_rows[c["name"]] = [c["name"]] + ([round(v * 1000, 4) for v in loc] if has else ["n/a"] * 3) + ["n/a", c["group"] or "n/a", c["bst_type"]]
        del F, X, m
# electrodes (Brainstorm SCS, metres -> mm)
wtsv(os.path.join(D, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "brainstorm_type"], list(elec_rows.values()))
wjson(os.path.join(D, f"{sub}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "mm",
      "iEEGCoordinateSystemDescription": "Brainstorm subject coordinate system (SCS) of the release channel file: origin midway between LPA and RPA, X towards the nasion, Y towards LPA, Z up; defined on the subject's own T1 MRI (fiducials in the release MRI file). Native space, converted from metres to millimetres; contacts typed SEEG_NO_LOC/EEG/ECG/REF have no location (n/a). The MRI is not redistributed (not defaced); the Zenodo record contains it (anat/Subject01/subjectimage_MRI.mat) together with cortex surfaces and volume atlases.",
      "iEEGCoordinateProcessingDescription": "as released (authors' electrode localisation in Brainstorm)", "iEEGCoordinateProcessingReference": "https://doi.org/10.1002/acn3.70115"})
wjson(os.path.join(OUT, "space-Other_electrodes.json"), {"SpatialReference": "n/a: the subject's T1 MRI (release file anat/Subject01/subjectimage_MRI.mat, 256x256x256, 1 mm) is not redistributed because it is not defaced",
      "brainstorm_type": {"Description": "Brainstorm channel type (SEEG, SEEG_NO_LOC = sEEG contact without location, EEG = scalp, EEG REF, ECG)"}})
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "group", "diagnosis"], [[sub, 25, "F", "epilepsy", "bilateral perisylvian polymicrogyria, intractable focal seizures (Zenodo description)"]])
wjson(os.path.join(OUT, "participants.json"), {"age": {"Description": "age in years (Zenodo description: 'A 25-year-old female')", "Units": "year"}, "sex": {"Description": "F = female (Zenodo description)"},
      "group": {"Description": "patient"}, "diagnosis": {"Description": "from the Zenodo description"}})
wtsv(os.path.join(OUT, sub, f"{sub}_scans.tsv"), ["filename", "acq_time"], [[f"ieeg/{r['stem']}_ieeg.vhdr", "n/a"] for r in rep["runs"]])
# sourcedata: Brainstorm data/channel/study files of the converted segments (no anat, no MRI), the protocol file and the two docx
for folder, fn, bfolder, bfn, task, tdesc in SEG:
    for fo in (folder, bfolder):
        if fo is None:
            continue
        for p in sorted(glob.glob(os.path.join(B, "data", "Subject01", fo, "*.mat"))):
            if os.path.basename(p).startswith(("data_", "channel", "brainstormstudy")):
                c = copy_verbatim(p, os.path.join(SD, "data", "Subject01", fo, os.path.basename(p))); rep["sourcedata"].append({"dst": os.path.relpath(c["dst"], OUT), "bytes": c["bytes"], "sha256": c["sha256"]})
for p in [os.path.join(B, "data", "protocol.mat")] + glob.glob(os.path.join(SRC, "files", "*.docx")):
    if os.path.exists(p):
        c = copy_verbatim(p, os.path.join(SD, os.path.basename(p))); rep["sourcedata"].append({"dst": os.path.relpath(c["dst"], OUT), "bytes": c["bytes"], "sha256": c["sha256"]})
wjson_any(os.path.join(SD, "intra_electrodes.json"), {"IntraElectrodes": ie_all, "note": "Brainstorm IntraElectrodes definitions of the release channel file (names, contact number/spacing/diameter/length, electrode length; metres)"})
hours = sum(r["duration_s"] for r in rep["runs"] if r["acq"] == "monopolar") / 3600
wjson(os.path.join(OUT, "dataset_description.json"), dataset_description(
    "VakilnaSEEG-fingerprint: sEEG seizure-onset and interictal segments of one polymicrogyria patient with electrode coordinates (Vakilna et al., 2025) (derivative)",
    authors, "CC-BY-4.0", DOI, URL, [f"https://doi.org/{DOI}", "https://doi.org/10.1002/acn3.70115"],
    f"Cite Vakilna YS et al. (2025) Time-Frequency Fingerprint Analysis in SEEG Source-Space to Identify the Epileptogenic Zone, Ann Clin Transl Neurol, doi:10.1002/acn3.70115, and the Zenodo dataset doi:{DOI}.",
    ["Collected under the source study's approval (University of Texas Health Science Center at Houston); see the article"],
    [{"Name": "o5c_B_conv_lit187.py", "Description": "Brainstorm sensor-space data files (notch-filtered, resampled to 1000 Hz by the authors) repackaged as BrainVision float32 in microvolts; Brainstorm events to events.tsv; channel locations (SCS) to electrodes.tsv. No further processing. Code in code/."}],
    dataset_type="derivative", keywords=list(meta.get("keywords", [])) + ["stereo-EEG", "epileptogenic zone", "seizure onset", "Brainstorm"],
    source_version=f"Zenodo version {meta.get('version', 'n/a')} ({meta.get('publication_date', 'n/a')})"))
desc = re.sub("<[^>]+>", " ", meta.get("description", "")).strip()
lines = "\n".join(f"  - `{r['stem']}`: {r['source']} ({r['comment']}), {r['n_channels']} ch, {r['sfreq']:g} Hz, {r['duration_s']:.0f} s, time axis {r['time_axis'][0]:.0f}..{r['time_axis'][1]:.0f} s, {r['n_events']} events" for r in rep["runs"])
open(os.path.join(OUT, "README"), "w").write(f"""# sEEG seizure-onset and interictal segments of one polymicrogyria patient, with electrode coordinates (Vakilna et al., 2025) - derivative

Release: Zenodo {ZEN} "{meta.get('title', '')}" ({', '.join(authors)}; doi:{DOI}, CC-BY-4.0), the Brainstorm database
(`BST_dataset.zip`) behind the case study doi:10.1002/acn3.70115. Zenodo description: "{desc}"

## Content (1 participant, {len(rep['runs'])} recordings, {hours * 60:.1f} min of sensor-space data)
The release contains three short sensor-space segments of the clinical stereo-EEG (Nihon Kohden EEG-1200A, 2000 Hz) that
the authors imported into Brainstorm, notch-filtered (180-900 Hz harmonics) and resampled to 1000 Hz, and their bipolar
('bipolar 2') montage versions. They are written here as BrainVision float32 in microvolts (Brainstorm stores volts),
`acq-monopolar` (267 channels: 228 located sEEG contacts, 20 sEEG contacts without location, scalp EEG, ECG, DC and
marker channels, as typed in the Brainstorm channel file) and `acq-bipolar` (206 bipolar pairs). This is a DERIVATIVE
dataset; the full raw recording is not in the release. Brainstorm events (clinical marks, the authors' onset/wave/spike
markers and channel annotations) are in events.tsv; the time axis of each segment (0 = reference event) is kept in the
`brainstorm_time` column.
{lines}

## Localisation
`electrodes.tsv` carries the Brainstorm channel locations (subject coordinate system SCS, native space, mm) of 228 sEEG
contacts of 16 depth electrodes ({', '.join(x['Name'] for x in (ie_all or []))}); electrode geometry (contact number,
spacing, diameter) is in `sourcedata/zenodo-{ZEN}/intra_electrodes.json`. The subject's T1 MRI, head/cortex surfaces and
volume atlases of the release are NOT redistributed (the MRI is not defaced); they, and the authors' source-space results
(OpenMEEG head model, sLORETA kernels, time-frequency maps, scout time series), remain in the Zenodo record.

## Provenance
Converted by the ieeg-nemar-20261005 campaign (lane o5c-B) with `code/o5c_B_conv_lit187.py`; round-trip checks (MNE
read-back) in code/conversion_report.json. sourcedata holds the converted Brainstorm data/channel/study files verbatim.
""")
copy_code(CID, ["o5c_B_conv_lit187.py", "o5c_B_common.py", "b3w10_common.py", "b3w3_common.py"])
bidsignore(OUT, ["sourcedata/", "code/"])
rep["hours"] = hours; rep["subjects"] = 1; rep["n_recordings"] = len(rep["runs"]); rep["n_located_contacts"] = sum(1 for v in elec_rows.values() if v[1] != "n/a")
rep["roundtrip_all_ok"] = all(r["roundtrip_ok"] for r in rep["runs"])
report(CID, rep)
print(json.dumps({"runs": len(rep["runs"]), "hours": hours, "rt": rep["roundtrip_all_ok"], "located": rep["n_located_contacts"]}, indent=1, default=str))
print("CONVERT_OK" if rep["roundtrip_all_ok"] else "CONVERT_ROUNDTRIP_FAIL")
