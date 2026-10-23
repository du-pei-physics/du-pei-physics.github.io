# Copied from claude/paper/fc_modular_check.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Numerical check of the modular step NS -> R used in field_content.md 6.4:
t(1 - 1/tau) = -(theta2/eta)^12(tau), and sum_j a_j t(1-1/tau)^j = Z^R(tau) (q-series from fc_closed_form)."""
import json
import mpmath as mp
import fc_closed_form as F
mp.mp.dps = 60
eta = lambda tau: mp.exp(2j * mp.pi * tau / 24) * mp.qp(mp.exp(2j * mp.pi * tau))
th = lambda k, tau: mp.jtheta(k, 0, mp.exp(1j * mp.pi * tau))
t = lambda tau: (th(3, tau) / eta(tau)) ** 12
for tau in (mp.mpc(0, 1.1), mp.mpc(0.3, 0.9), mp.mpc(-0.2, 1.5)):
    print("t(1-1/tau) + (th2/eta)^12, rel. err:", mp.nstr(abs(t(1 - 1 / tau) + (th(2, tau) / eta(tau)) ** 12) / abs(t(1 - 1 / tau)), 3))
d = json.load(open(F.HERE / "out" / "fc_closed_form.json"))
for name in ("T", "Tprime", "P"):
    a = d[name]["a_j"]; ZR = F.eval_poly(a, "R")
    for tau in (mp.mpc(0.1, 1.6), mp.mpc(-0.3, 2.2)):
        q = mp.exp(2j * mp.pi * tau)
        mod = sum(c * t(1 - 1 / tau) ** j for j, c in enumerate(a))
        ser = sum(v * q ** (k // 2) for k, v in ZR.items() if k <= 50)
        print(name, tau, "rel. err:", mp.nstr(abs(mod - ser) / abs(ser), 3))
