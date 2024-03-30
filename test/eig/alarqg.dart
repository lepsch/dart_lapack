import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

Future<void> alarqg(
  final String PATH,
  final int NMATS,
  final Array<bool> DOTYPE_,
  final int NTYPES,
  final Nin NIN,
  final Nout NOUT,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  bool FIRSTT;
  String C1;
  String LINE;
  int I, I1, IC = 0, J, K, LENP, NT;
  final NREQ = Array<int>(100);
  const INTSTR = '0123456789';

  if (NMATS >= NTYPES) {
    // Test everything if NMATS >= NTYPES.

    for (I = 1; I <= NTYPES; I++) {
      DOTYPE[I] = true;
    }
  } else {
    for (I = 1; I <= NTYPES; I++) {
      DOTYPE[I] = false;
    }
    FIRSTT = true;

    // Read a line of matrix types if 0 < NMATS < NTYPES.
    try {
      if (NMATS > 0) {
        LINE = await NIN.readLine();
        LENP = LINE.length;
        I = 0;
        nextChar:
        for (J = 1; J <= NMATS; J++) {
          NREQ[J] = 0;
          I1 = 0;
          while (true) {
            I++;
            if (I > LENP) {
              if (J == NMATS && I1 > 0) {
                continue nextChar;
              } else {
                NOUT.println(
                    '\n\n *** Not enough matrix types on input line/n$LINE');
                print9994(NOUT, NMATS);
                return;
              }
            }
            if (LINE.substring(I - 1, I) != ' ' &&
                LINE.substring(I - 1, I) != ',') {
              I1 = I;
              C1 = LINE.substring(I1 - 1, I1);

              // Check that a valid integer was read
              var isDigit = false;
              for (K = 1; K <= 10; K++) {
                if (C1 == INTSTR.substring(K - 1, K)) {
                  IC = K - 1;
                  isDigit = true;
                  break;
                }
              }
              if (!isDigit) {
                NOUT.println(
                    '\n\n *** Invalid integer value in column ${I.i2} of input line:\n$LINE');
                print9994(NOUT, NMATS);
                return;
              }
              NREQ[J] = 10 * NREQ[J] + IC;
            } else if (I1 > 0) {
              continue nextChar;
            }
          }
        }
      }
      for (I = 1; I <= NMATS; I++) {
        NT = NREQ[I];
        if (NT > 0 && NT <= NTYPES) {
          if (DOTYPE[NT]) {
            if (FIRSTT) NOUT.println('');
            FIRSTT = false;
            NOUT.println(
                ' *** Warning:  duplicate request of matrix type ${NT.i2} for ${PATH.a3}');
          }
          DOTYPE[NT] = true;
        } else {
          NOUT.println(
              ' *** Invalid type request for $PATH, type  ${NT.i4}: must satisfy  1 <= type <= ${NTYPES.i2}');
        }
      }
    } on EOF catch (_) {
      NOUT.println(
          '\n *** End of file reached when trying to read matrix types for ${PATH.a3}\n *** Check that you are requesting the right number of types for each path\n');

      NOUT.println();
      rethrow;
    }
  }
}

void print9994(final Nout NOUT, final int nmats) {
  NOUT.println(
      ' ==> Specify ${nmats.i4} matrix types on this line or adjust NTYPES on previous line');
}
