#pragma once

// Read an EMTSim netlist. One element per line; node 0 = ground; '#' starts a comment,
// and the first word of the comment (if any) names the element in the recorded currents.
// Optional parameters can be given as key=value. The same format is read by the Python
// version (netlist.py).
//
//   R  a b  ohm                     L  a b  henry                 C  a b  farad
//   V  a b  peak angle freq [R] [L]     R > 0: Norton; L: ideal source behind series R + L
//   SW a b  closed|open t1 t2 ...  [czero] [chop=A] [ron=] [roff=]
//   F  a b  t_apply [t_clear]  [ron=] [roff=]
//   RL3 a1 b1 c1 a2 b2 c2  R1 L1 [R0 L0]            coupled three-phase series R-L (sequence data)
//   PI  a b  R L C                                   single-phase nominal-pi line (C total)
//   PI3 a1 b1 c1 a2 b2 c2  R1 L1 C1 [R0 L0 C0]       three-phase nominal-pi line
//   V3  a b c  peak angle freq [R1 L1 [R0 L0]]       three-phase source (behind impedance if given)
//   TL  k m  Zc tau [R]                              Bergeron line
//   TL3 ka kb kc ma mb mc  Z1 tau1 Z0 tau0 [R1 R0]   three-phase Bergeron line (transposed)
//   XF  p q s t  S V1 V2 f x r  [im=%] [knee=pu lsat=pu]           single-phase transformer
//   XF3 a1 b1 c1 a2 b2 c2  conn S V1 V2 f x r  [im=%] [knee= lsat=]  conn: YNyn0 YNd1 YNd11 Dyn1 Dyn11 Yy0 Dd0 ...
//   LSAT a b  i1 flux1 i2 flux2 ...                  saturable inductor
//   MOV a b  Vref [Iref] [q]                         surge arrester
//   D   a k  [ron=] [roff=] [vf=]                    diode
//   THY a k  alpha freq [ref=deg] [width=deg] [ron=] [roff=]   thyristor
//   SM  a b c  S V f H xd2 ra E delta0 Pm [D=]       synchronous machine (classical)
//   IM  a b c  S V f poles rs xls rr xlr xm H Tload [speed0=]   induction machine
//   .tran dt tstop                                   time step and stop time (before the elements)

#include <string>

#include "circuit.h"

Circuit loadNetlist(const std::string& path);   // throws std::runtime_error with file:line on errors
Circuit parseNetlist(const std::string& text, const std::string& where = "netlist");
int parseOnOff(const std::string& value);   // on/1/true/yes -> 1, off/0/false/no -> 0, otherwise -1
