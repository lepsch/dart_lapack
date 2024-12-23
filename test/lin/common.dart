// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:async/async.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

class _Cenvir {
  int NPROC = 0;
  int NSHIFT = 0;
  int MAXB = 0;
}

final cenvir = _Cenvir();

class _Claenv {
  Array<int> IPARMS = Array<int>(100);
}

final claenv = _Claenv();

class _Sslct {
  int SELOPT = 0;
  int SELDIM = 0;
  Array<bool> SELVAL = Array<bool>(20);
  Array<double> SELWR = Array<double>(20);
  Array<double> SELWI = Array<double>(20);
}

final sslct = _Sslct();

class _Mn {
  int M = 0;
  int N = 0;
  int MPLUSN = 0;
  int K = 0;
  bool FS = false;
}

final mn = _Mn();

class _Srnamc {
  String SRNAMT = '';
}

final srnamc = _Srnamc();

class _Infoc {
  int INFOT = 0;
  Nout NOUT = Nout(NullStreamSink());
  Box<bool> OK = Box(false);
  Box<bool> LERR = Box(false);
}

final infoc = _Infoc();
