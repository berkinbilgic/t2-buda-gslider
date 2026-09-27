# T2-BUDA-gSlider: Efficient T2 Mapping with Blip-Up/Down EPI and gSlider-SMS

MATLAB code and example in vivo data for rapid high-resolution T2 mapping with blip-up/down acquisition (BUDA), gSlider-SMS and subspace reconstruction.

`github_001_lite_release.m` reconstructs the blip-up/down EPI data with a MUSSELS-type joint reconstruction that uses a B0 field map model, performs the gSlider slice reconstruction with the simulated RF slice profiles, and estimates T2 maps with dictionary matching.

## Data download

`data_pack.mat` (178 MB) is larger than GitHub's 100 MB file limit, so it is attached to the [v1.0 release](https://github.com/berkinbilgic/t2-buda-gslider/releases/tag/v1.0) instead of the repository. Download it into the repository folder:

```bash
curl -LO https://github.com/berkinbilgic/t2-buda-gslider/releases/download/v1.0/data_pack.mat
```

## Usage

Run `github_001_lite_release.m` from the repository folder in MATLAB.

## Third-party code

`tools/` includes Michael Lustig's SparseMRI and coil compression code (Zhang et al., MRM 2013; see `tools/README`), MRI reconstruction utilities by Jonathan Polimeni (`mrir_*`), and WaveLab wavelet MEX files (copyright notices retained in the files).

## Reference

X Cao, C Liao, Z Zhang, SS Iyer, K Wang, H He, H Liu, K Setsompop, J Zhong, B Bilgic. Efficient T2 mapping with Blip-up/down EPI and gSlider-SMS (T2-BUDA-gSlider). [arXiv:1909.12999](https://arxiv.org/abs/1909.12999)

Contact: berkin AT nmr.mgh.harvard.edu
