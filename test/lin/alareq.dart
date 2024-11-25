// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

Future<void> alareq(
  final String PATH,
  final int NMATS,
  final Array<bool> DOTYPE_,
  final int NTYPES,
  final Nin NIN,
  final Nout NOUT,
) async {
  final DOTYPE = DOTYPE_.having();
  final NREQ = Array<int>(100);
  const INTSTR = '0123456789';

  if (NMATS >= NTYPES) {
    // Test everything if NMATS >= NTYPES.

    for (var I = 1; I <= NTYPES; I++) {
      DOTYPE[I] = true;
    }
  } else {
    for (var I = 1; I <= NTYPES; I++) {
      DOTYPE[I] = false;
    }

    // Read a line of matrix types if 0 < NMATS < NTYPES.

    if (NMATS > 0) {
      final String LINE;
      try {
        LINE = await NIN.readLine();
      } on EOF catch (_) {
        NOUT.println(
            '\n *** End of file reached when trying to read matrix types for ${PATH.a3}\n *** Check that you are requesting the right number of types for each path\n');
        NOUT.println();
        rethrow;
      }
      final LENP = LINE.length;
      var I = 0, IC = 0;
      nextValue:
      for (var J = 1; J <= NMATS; J++) {
        NREQ[J] = 0;
        var I1 = 0;
        nextChar:
        while (true) {
          I++;
          if (I > LENP) {
            if (J == NMATS && I1 > 0) {
              continue nextValue;
            } else {
              NOUT.println(
                  '\n\n *** Not enough matrix types on input line\n${LINE.a79}');
              NOUT.printInfo(NMATS);
              return;
            }
          }
          if (LINE[I - 1] != ' ' && LINE[I - 1] != ',') {
            I1 = I;
            final C1 = LINE[I1 - 1];

            // Check that a valid integer was read
            var isValidInt = false;
            for (var K = 1; K <= 10; K++) {
              if (C1 == INTSTR[K - 1]) {
                IC = K - 1;
                isValidInt = true;
                break;
              }
            }

            if (!isValidInt) {
              NOUT.println(
                  '\n\n *** Invalid integer value in column ${I.i2} of input line:\n${LINE.a79}');
              NOUT.printInfo(NMATS);
              return;
            }
            NREQ[J] = 10 * NREQ[J] + IC;
            continue nextChar;
          } else if (I1 > 0) {
            continue nextValue;
          }
        }
      }
    }
    var FIRSTT = true;
    for (var I = 1; I <= NMATS; I++) {
      final NT = NREQ[I];
      if (NT > 0 && NT <= NTYPES) {
        if (DOTYPE[NT]) {
          if (FIRSTT) NOUT.println();
          FIRSTT = false;
          NOUT.println(
              ' *** Warning:  duplicate request of matrix type ${NT.i2} for ${PATH.a3}');
        }
        DOTYPE[NT] = true;
      } else {
        NOUT.println(
            ' *** Invalid type request for ${PATH.a3}, type  ${NT.i4}: must satisfy  1 <= type <= ${NTYPES.i2}');
      }
    }
  }
}

extension on Nout {
  void printInfo(int NMATS) {
    println(
        ' ==> Specify ${NMATS.i4} matrix types on this line or adjust NTYPES on previous line');
  }
}
