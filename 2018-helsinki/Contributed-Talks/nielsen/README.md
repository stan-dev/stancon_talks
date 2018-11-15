# Paper for submission to StanCon

## Build instructions

1. clone this repository (stancon)
2. Install Stan (see below)
3. Install Haskell Stack (https://docs.haskellstack.org/en/stable/install_and_upgrade/)
4. `cd stancon`
5. `stack install`
6. `stancon >stancon.html`
7. now you can open `stancon.html` in a web browser. For instance, you could enter at the command line: `firefox stancon.html’

The paper is using the inliterate system (https://github.com/diffusionkinetics/open/tree/master/inliterate)

## Installing Stan

you will need various compiler infrastructure, for instance clang and libc++-dev.

```
wget https://github.com/stan-dev/cmdstan/releases/download/v2.17.0/cmdstan-2.17.0.tar.gz
tar -xzvf cmdstan-2.17.0.tar.gz
sudo mv cmdstan-2.17.0 /opt/stan
(cd /opt/stan && sudo make build -j4)
rm cmdstan-2.17.0.tar.gz
```

If you have installed Stan in a different directory than /opt/stan, then set the environment variable CMDSTAN_HOME must point to the directory in which you have installed cmdstan
