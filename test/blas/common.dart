// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:async/async.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/nio.dart';

class _Combla {
  int ICASE = 0;
  int N = 0;
  int INCX = 0;
  int INCY = 0;
  int MODE = 0;
  bool PASS = false;
}

final combla = _Combla();

class _Infoc {
  int INFOT = 0;
  Nout NOUT = Nout(NullStreamSink());
  Box<bool> OK = Box(false);
  Box<bool> LERR = Box(false);

  Nout get NOUTC => NOUT;
  set NOUTC(Nout noutc) => NOUT = noutc;
}

final infoc = _Infoc();

class _Srnamc {
  String SRNAMT = '';
}

final srnamc = _Srnamc();
