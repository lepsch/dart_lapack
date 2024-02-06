import 'dart:async';

import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

Future<void> alareq(
  final String PATH,
  final int NMATS,
  final Array<bool> DOTYPE,
  final int NTYPES,
  final Nin NIN,
  final Nout NOUT,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool FIRSTT;
  String C1;
  String LINE = '';
  int I, I1, IC = 0, J, K, LENP, NT;
  final NREQ = Array<int>(100);
  const INTSTR = '0123456789';

  try {
    if (NMATS >= NTYPES) {
      // Test everything if NMATS >= NTYPES.

      for (I = 1; I <= NTYPES; I++) {
        DOTYPE[I] = true;
      }
      return;
    }

    for (I = 1; I <= NTYPES; I++) {
      DOTYPE[I] = false;
    }
    FIRSTT = true;

    // Read a line of matrix types if 0 < NMATS < NTYPES.

    if (NMATS > 0) {
      LINE = await NIN.readLine();
      LENP = LINE.length;
      I = 0;
      nextValue:
      for (J = 1; J <= NMATS; J++) {
        NREQ[J] = 0;
        I1 = 0;
        nextChar:
        while (true) {
          I = I + 1;
          if (I > LENP) {
            if (J == NMATS && I1 > 0) {
              continue nextValue;
            }

            print9995(NOUT, LINE);
            print9994(NOUT, NMATS);
            return;
          }

          if (LINE.substring(I - 1, I) != ' ' &&
              LINE.substring(I - 1, I) != ',') {
            I1 = I;
            C1 = LINE.substring(I1 - 1, I1);

            // Check that a valid integer was read
            var isValidInt = false;
            for (K = 1; K <= 10; K++) {
              if (C1 == INTSTR.substring(K - 1, K)) {
                IC = K - 1;
                isValidInt = true;
                break;
              }
            }
            if (!isValidInt) {
              print9996(NOUT, I, LINE);
              print9994(NOUT, NMATS);
              return;
            }

            NREQ[J] = 10 * NREQ[J] + IC;
            continue nextChar;
          } else if (I1 > 0) {
            continue;
          } else {
            continue nextChar;
          }
        }
      }
    }
    for (I = 1; I <= NMATS; I++) {
      NT = NREQ[I];
      if (NT > 0 && NT <= NTYPES) {
        if (DOTYPE[NT]) {
          if (FIRSTT) NOUT.println();
          FIRSTT = false;
          print9997(NOUT, NT, PATH);
        }
        DOTYPE[NT] = true;
      } else {
        print9999(NOUT, PATH, NT, NTYPES);
      }
    }
  } on EOF catch (_) {
    NOUT.println(
      '\n *** End of file reached when trying to read matrix types for ${PATH.a3}\n *** Check that you are requesting the right number of types for each path\n',
    );
    NOUT.println();
    rethrow;
  }
}

void print9999(
  final Nout nout,
  final String path,
  final int nt,
  final int ntypes,
) {
  nout.println(
    ' *** Invalid type request for ${path.a3}, type  ${nt.i4}: must satisfy  1 <= type <= ${ntypes.i2}',
  );
}

void print9997(final Nout nout, final int nt, final String path) {
  nout.println(
    ' *** Warning:  duplicate request of matrix type ${nt.i2} for ${path.a3}',
  );
}

void print9996(final Nout nout, final int column, final String line) {
  nout.println(
    ' *** Invalid integer value in column ${column.i2} of input line:\n ${line.a79}',
  );
}

void print9995(final Nout nout, final String line) {
  nout.println(' *** Not enough matrix types on input line\n${line.a79}');
}

void print9994(final Nout nout, final int n) {
  nout.println(
    ' ==> Specify ${n.i4} matrix types on this line or adjust NTYPES on previous line',
  );
}
