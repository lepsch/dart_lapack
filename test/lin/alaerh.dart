// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/nio.dart';

import 'aladhd.dart';
import 'alahd.dart';

void alaerh(
  final String PATH,
  final String SUBNAM,
  final int INFO,
  final int INFOE,
  final String OPTS,
  final int M,
  final int N,
  final int KL,
  final int KU,
  final int N5,
  final int IMAT,
  final int NFAIL,
  final Box<int> NERRS,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  String UPLO;

  if (INFO == 0) return;
  final P2 = PATH.substring(1, 3);
  final C3 = SUBNAM.substring(3, 6);

  // Print the header if this is the first error message.

  if (NFAIL == 0 && NERRS.value == 0) {
    if (lsamen(3, C3, 'SV ') || lsamen(3, C3, 'SVX')) {
      aladhd(NOUT, PATH);
    } else {
      alahd(NOUT, PATH);
    }
  }
  NERRS.value++;

  // Print the message detailing the error and form of recovery,
  // if any.

  if (lsamen(2, P2, 'GE')) {
    // xGE:  General matrices

    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9988(SUBNAM.trim(), INFO, INFOE, M, N, N5, IMAT);
      } else {
        NOUT.print9975(SUBNAM.trim(), INFO, M, N, N5, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9984(SUBNAM.trim(), INFO, INFOE, N, N5, IMAT);
      } else {
        NOUT.print9970(SUBNAM.trim(), INFO, N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9992(
            SUBNAM.trim(), INFO, INFOE, OPTS[0], OPTS[1], N, N5, IMAT);
      } else {
        NOUT.print9997(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'TRI')) {
      NOUT.print9971(SUBNAM.trim(), INFO, N, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS')) {
      NOUT.print9978(SUBNAM.trim(), INFO, M, N, IMAT);
    } else if (lsamen(3, C3, 'CON')) {
      NOUT.print9969(SUBNAM.trim(), INFO, OPTS[0], M, IMAT);
    } else if (lsamen(3, C3, 'LS ')) {
      NOUT.print9965(SUBNAM.trim(), INFO, OPTS[0], M, N, KL, N5, IMAT);
    } else if (lsamen(3, C3, 'LSX') || lsamen(3, C3, 'LSS')) {
      NOUT.print9974(SUBNAM.trim(), INFO, M, N, KL, N5, IMAT);
    } else {
      NOUT.print9963(SUBNAM.trim(), INFO, OPTS[0], M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'GB')) {
    // xGB:  General band matrices

    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9989(SUBNAM.trim(), INFO, INFOE, M, N, KL, KU, N5, IMAT);
      } else {
        NOUT.print9976(SUBNAM.trim(), INFO, M, N, KL, KU, N5, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9986(SUBNAM.trim(), INFO, INFOE, N, KL, KU, N5, IMAT);
      } else {
        NOUT.print9972(SUBNAM.trim(), INFO, N, KL, KU, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9993(
            SUBNAM.trim(), INFO, INFOE, OPTS[0], OPTS[1], N, KL, KU, N5, IMAT);
      } else {
        NOUT.print9998(
            SUBNAM.trim(), INFO, OPTS[0], OPTS[1], N, KL, KU, N5, IMAT);
      }
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS')) {
      NOUT.print9977(SUBNAM.trim(), INFO, M, N, KL, KU, IMAT);
    } else if (lsamen(3, C3, 'CON')) {
      NOUT.print9968(SUBNAM.trim(), INFO, OPTS[0], M, KL, KU, IMAT);
    } else {
      NOUT.print9964(SUBNAM.trim(), INFO, OPTS[0], M, KL, KU, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'GT')) {
    // xGT:  General tridiagonal matrices

    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9987(SUBNAM.trim(), INFO, INFOE, N, IMAT);
      } else {
        NOUT.print9973(SUBNAM.trim(), INFO, N, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9984(SUBNAM.trim(), INFO, INFOE, N, N5, IMAT);
      } else {
        NOUT.print9970(SUBNAM.trim(), INFO, N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9992(
            SUBNAM.trim(), INFO, INFOE, OPTS[0], OPTS[1], N, N5, IMAT);
      } else {
        NOUT.print9997(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'CON')) {
      NOUT.print9969(SUBNAM.trim(), INFO, OPTS[0], M, IMAT);
    } else {
      NOUT.print9963(SUBNAM.trim(), INFO, OPTS[0], M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'PO')) {
    // xPO:  Symmetric or Hermitian positive definite matrices

    UPLO = OPTS[0];
    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9980(SUBNAM.trim(), INFO, INFOE, UPLO, M, N5, IMAT);
      } else {
        NOUT.print9956(SUBNAM.trim(), INFO, UPLO, M, N5, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9979(SUBNAM.trim(), INFO, INFOE, UPLO, N, N5, IMAT);
      } else {
        NOUT.print9955(SUBNAM.trim(), INFO, UPLO, N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9990(
            SUBNAM.trim(), INFO, INFOE, OPTS[0], OPTS[1], N, N5, IMAT);
      } else {
        NOUT.print9995(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'TRI')) {
      NOUT.print9956(SUBNAM.trim(), INFO, UPLO, M, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS') ||
        lsamen(3, C3, 'CON')) {
      NOUT.print9960(SUBNAM.trim(), INFO, UPLO, M, IMAT);
    } else {
      NOUT.print9955(SUBNAM.trim(), INFO, UPLO, M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'PS')) {
    // xPS:  Symmetric or Hermitian positive semi-definite matrices

    UPLO = OPTS[0];
    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9980(SUBNAM, INFO, INFOE, UPLO, M, N5, IMAT);
      } else {
        NOUT.print9956(SUBNAM, INFO, UPLO, M, N5, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9979(SUBNAM, INFO, INFOE, UPLO, N, N5, IMAT);
      } else {
        NOUT.print9955(SUBNAM, INFO, UPLO, N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9990(SUBNAM, INFO, INFOE, OPTS[0], OPTS[1], N, N5, IMAT);
      } else {
        NOUT.print9995(SUBNAM, INFO, OPTS[0], OPTS[1], N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'TRI')) {
      NOUT.print9956(SUBNAM, INFO, UPLO, M, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMT') ||
        lsamen(3, C3, 'CON')) {
      NOUT.print9960(SUBNAM, INFO, UPLO, M, IMAT);
    } else {
      NOUT.print9955(SUBNAM, INFO, UPLO, M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'SY') ||
      lsamen(2, P2, 'SR') ||
      lsamen(2, P2, 'SK') ||
      lsamen(2, P2, 'HE') ||
      lsamen(2, P2, 'HR') ||
      lsamen(2, P2, 'HK') ||
      lsamen(2, P2, 'HA')) {
    // xSY: symmetric indefinite matrices
    //      with partial (Bunch-Kaufman) pivoting;
    // xSR: symmetric indefinite matrices
    //      with rook (bounded Bunch-Kaufman) pivoting;
    // xSK: symmetric indefinite matrices
    //      with rook (bounded Bunch-Kaufman) pivoting,
    //      new storage format;
    // xHE: Hermitian indefinite matrices
    //      with partial (Bunch-Kaufman) pivoting.
    // xHR: Hermitian indefinite matrices
    //      with rook (bounded Bunch-Kaufman) pivoting;
    // xHK: Hermitian indefinite matrices
    //      with rook (bounded Bunch-Kaufman) pivoting,
    //      new storage format;
    // xHA: Hermitian matrices
    //      Aasen Algorithm

    UPLO = OPTS[0];
    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9980(SUBNAM.trim(), INFO, INFOE, UPLO, M, N5, IMAT);
      } else {
        NOUT.print9956(SUBNAM.trim(), INFO, UPLO, M, N5, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(2, C3, 'SV')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9979(SUBNAM.trim(), INFO, INFOE, UPLO, N, N5, IMAT);
      } else {
        NOUT.print9955(SUBNAM.trim(), INFO, UPLO, N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9990(
            SUBNAM.trim(), INFO, INFOE, OPTS[0], OPTS[1], N, N5, IMAT);
      } else {
        NOUT.print9995(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], N, N5, IMAT);
      }
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS') ||
        lsamen(3, C3, 'TRI') ||
        lsamen(3, C3, 'CON')) {
      NOUT.print9960(SUBNAM.trim(), INFO, UPLO, M, IMAT);
    } else {
      NOUT.print9955(SUBNAM.trim(), INFO, UPLO, M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'PP') ||
      lsamen(2, P2, 'SP') ||
      lsamen(2, P2, 'HP')) {
    // xPP, xHP, or xSP:  Symmetric or Hermitian packed matrices

    UPLO = OPTS[0];
    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9983(SUBNAM.trim(), INFO, INFOE, UPLO, M, IMAT);
      } else {
        NOUT.print9960(SUBNAM.trim(), INFO, UPLO, M, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9979(SUBNAM.trim(), INFO, INFOE, UPLO, N, N5, IMAT);
      } else {
        NOUT.print9955(SUBNAM.trim(), INFO, UPLO, N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9990(
            SUBNAM.trim(), INFO, INFOE, OPTS[0], OPTS[1], N, N5, IMAT);
      } else {
        NOUT.print9995(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], N, N5, IMAT);
      }
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS') ||
        lsamen(3, C3, 'TRI') ||
        lsamen(3, C3, 'CON')) {
      NOUT.print9960(SUBNAM.trim(), INFO, UPLO, M, IMAT);
    } else {
      NOUT.print9955(SUBNAM.trim(), INFO, UPLO, M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'PB')) {
    // xPB:  Symmetric (Hermitian) positive definite band matrix

    UPLO = OPTS[0];
    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9982(SUBNAM.trim(), INFO, INFOE, UPLO, M, KL, N5, IMAT);
      } else {
        NOUT.print9958(SUBNAM.trim(), INFO, UPLO, M, KL, N5, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9981(SUBNAM.trim(), INFO, INFOE, UPLO, N, KL, N5, IMAT);
      } else {
        NOUT.print9957(SUBNAM.trim(), INFO, UPLO, N, KL, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9991(
            SUBNAM.trim(), INFO, INFOE, OPTS[0], OPTS[1], N, KL, N5, IMAT);
      } else {
        NOUT.print9996(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], N, KL, N5, IMAT);
      }
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS') ||
        lsamen(3, C3, 'CON')) {
      NOUT.print9959(SUBNAM.trim(), INFO, UPLO, M, KL, IMAT);
    } else {
      NOUT.print9957(SUBNAM.trim(), INFO, UPLO, M, KL, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'PT')) {
    // xPT:  Positive definite tridiagonal matrices

    if (lsamen(3, C3, 'TRF')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9987(SUBNAM.trim(), INFO, INFOE, N, IMAT);
      } else {
        NOUT.print9973(SUBNAM.trim(), INFO, N, IMAT);
      }
      if (INFO != 0) NOUT.print9949();
    } else if (lsamen(3, C3, 'SV ')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9984(SUBNAM.trim(), INFO, INFOE, N, N5, IMAT);
      } else {
        NOUT.print9970(SUBNAM.trim(), INFO, N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'SVX')) {
      if (INFO != INFOE && INFOE != 0) {
        NOUT.print9994(SUBNAM.trim(), INFO, INFOE, OPTS[0], N, N5, IMAT);
      } else {
        NOUT.print9999(SUBNAM.trim(), INFO, OPTS[0], N, N5, IMAT);
      }
    } else if (lsamen(3, C3, 'CON')) {
      if (lsame(SUBNAM[0], 'S') || lsame(SUBNAM[0], 'D')) {
        NOUT.print9973(SUBNAM.trim(), INFO, M, IMAT);
      } else {
        NOUT.print9969(SUBNAM.trim(), INFO, OPTS[0], M, IMAT);
      }
    } else {
      NOUT.print9963(SUBNAM.trim(), INFO, OPTS[0], M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'TR')) {
    // xTR:  Triangular matrix

    if (lsamen(3, C3, 'TRI')) {
      NOUT.print9961(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], M, N5, IMAT);
    } else if (lsamen(3, C3, 'CON')) {
      NOUT.print9967(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], M, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATRS')) {
      NOUT.print9952(
          SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], OPTS[3], M, IMAT);
    } else {
      NOUT.print9953(
          SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'TP')) {
    // xTP:  Triangular packed matrix

    if (lsamen(3, C3, 'TRI')) {
      NOUT.print9962(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], M, IMAT);
    } else if (lsamen(3, C3, 'CON')) {
      NOUT.print9967(SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], M, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATPS')) {
      NOUT.print9952(
          SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], OPTS[3], M, IMAT);
    } else {
      NOUT.print9953(
          SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], M, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'TB')) {
    // xTB:  Triangular band matrix

    if (lsamen(3, C3, 'CON')) {
      NOUT.print9966(
          SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], M, KL, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATBS')) {
      NOUT.print9951(
          SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], OPTS[3], M, KL, IMAT);
    } else {
      NOUT.print9954(
          SUBNAM.trim(), INFO, OPTS[0], OPTS[1], OPTS[2], M, KL, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'QR')) {
    // xQR:  QR factorization

    if (lsamen(3, C3, 'QRS')) {
      NOUT.print9974(SUBNAM.trim(), INFO, M, N, KL, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS')) {
      NOUT.print9978(SUBNAM.trim(), INFO, M, N, IMAT);
    }
  } else if (lsamen(2, P2, 'QK')) {
    // xQK:  truncated QR factorization with pivoting

    if (lsamen(7, SUBNAM.substring(1, 8), 'GEQP3RK')) {
      NOUT.print9930(SUBNAM.trim(), INFO, M, N, KL, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS')) {
      NOUT.print9978(SUBNAM.trim(), INFO, M, N, IMAT);
    }
  } else if (lsamen(2, P2, 'LQ')) {
    // xLQ:  LQ factorization

    if (lsamen(3, C3, 'LQS')) {
      NOUT.print9974(SUBNAM.trim(), INFO, M, N, KL, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS')) {
      NOUT.print9978(SUBNAM.trim(), INFO, M, N, IMAT);
    }
  } else if (lsamen(2, P2, 'QL')) {
    // xQL:  QL factorization

    if (lsamen(3, C3, 'QLS')) {
      NOUT.print9974(SUBNAM.trim(), INFO, M, N, KL, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS')) {
      NOUT.print9978(SUBNAM.trim(), INFO, M, N, IMAT);
    }
  } else if (lsamen(2, P2, 'RQ')) {
    // xRQ:  RQ factorization

    if (lsamen(3, C3, 'RQS')) {
      NOUT.print9974(SUBNAM.trim(), INFO, M, N, KL, N5, IMAT);
    } else if (lsamen(5, SUBNAM.substring(1, 6), 'LATMS')) {
      NOUT.print9978(SUBNAM.trim(), INFO, M, N, IMAT);
    }
  } else if (lsamen(2, P2, 'LU')) {
    if (INFO != INFOE && INFOE != 0) {
      NOUT.print9988(SUBNAM.trim(), INFO, INFOE, M, N, N5, IMAT);
    } else {
      NOUT.print9975(SUBNAM.trim(), INFO, M, N, N5, IMAT);
    }
  } else if (lsamen(2, P2, 'CH')) {
    if (INFO != INFOE && INFOE != 0) {
      NOUT.print9985(SUBNAM.trim(), INFO, INFOE, M, N5, IMAT);
    } else {
      NOUT.print9971(SUBNAM.trim(), INFO, M, N5, IMAT);
    }
  } else {
    // Print a generic message if the path is unknown.

    NOUT.print9950(SUBNAM.trim(), INFO);
  }
}

extension on Nout {
  // Description of error message (alphabetical, left to right)

  void print9999(
      String SUBNAM, int INFO, String FACT, int N, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM=${INFO.i5}, FACT=\'${FACT.a1}\', N=${N.i5}, NRHS=${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9998(String SUBNAM, int INFO, String FACT, String TRANS, int N,
      int KL, int KU, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> FACT=\'${FACT.a1}\', TRANS=\'${TRANS.a1}\', N=${N.i5}, KL=${KL.i5}, KU=${KU.i5}, NRHS=${NRHS.i4}, type ${TYPE.i1}');
  }

  void print9997(String SUBNAM, int INFO, String FACT, String TRANS, int N,
      int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> FACT=\'${FACT.a1}\', TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9996(String SUBNAM, int INFO, String FACT, String UPLO, int N,
      int KD, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, KD=${KD.i5}, NRHS=${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9995(String SUBNAM, int INFO, String FACT, String UPLO, int N,
      int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9994(String SUBNAM, int INFO, int INFOE, String FACT, int N,
      int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> FACT=\'${FACT.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9993(String SUBNAM, int INFO, int INFOE, String FACT, String TRANS,
      int N, int KL, int KU, int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> FACT=\'${FACT.a1}\', TRANS=\'${TRANS.a1}\', N=${N.i5}, KL=${KL.i5}, KU=${KU.i5}, NRHS=${NRHS.i4}, type ${TYPE.i1}');
  }

  void print9992(String SUBNAM, int INFO, int INFOE, String FACT, String TRANS,
      int N, int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> FACT=\'${FACT.a1}\', TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9991(String SUBNAM, int INFO, int INFOE, String FACT, String UPLO,
      int N, int KD, int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, KD=${KD.i5}, NRHS=${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9990(String SUBNAM, int INFO, int INFOE, String FACT, String UPLO,
      int N, int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9989(String SUBNAM, int INFO, int INFOE, int M, int N, int KL,
      int KU, int NB, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> M = ${M.i5}, N =${N.i5}, KL =${KL.i5}, KU =${KU.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9988(
      String SUBNAM, int INFO, int INFOE, int M, int N, int NB, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> M =${M.i5}, N =${N.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9987(String SUBNAM, int INFO, int INFOE, int N, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2} for N=${N.i5}, type ${TYPE.i2}');
  }

  void print9986(String SUBNAM, int INFO, int INFOE, int N, int KL, int KU,
      int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> N =${N.i5}, KL =${KL.i5}, KU =${KU.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9985(String SUBNAM, int INFO, int INFOE, int N, int NB, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> N =${N.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9984(
      String SUBNAM, int INFO, int INFOE, int N, int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9983(
      String SUBNAM, int INFO, int INFOE, String UPLO, int N, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, type ${TYPE.i2}');
  }

  void print9982(String SUBNAM, int INFO, int INFOE, String UPLO, int N, int KD,
      int NB, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, KD =${KD.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9981(String SUBNAM, int INFO, int INFOE, String UPLO, int N, int KD,
      int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> UPLO=\'${UPLO.a1}\', N =${N.i5}, KD =${KD.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9980(String SUBNAM, int INFO, int INFOE, String UPLO, int N, int NB,
      int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9979(String SUBNAM, int INFO, int INFOE, String UPLO, int N,
      int NRHS, int TYPE) {
    println(
        ' *** $SUBNAM returned with INFO =${INFO.i5} instead of ${INFOE.i2}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9978(String SUBNAM, int INFO, int M, int N, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5} for M =${M.i5}, N =${N.i5}, type ${TYPE.i2}');
  }

  void print9977(
      String SUBNAM, int INFO, int M, int N, int KL, int KU, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> M = ${M.i5}, N =${N.i5}, KL =${KL.i5}, KU =${KU.i5}, type ${TYPE.i2}');
  }

  void print9976(
      String SUBNAM, int INFO, int M, int N, int KL, int KU, int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> M = ${M.i5}, N =${N.i5}, KL =${KL.i5}, KU =${KU.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9975(String SUBNAM, int INFO, int M, int N, int NB, int IMAT) {
    println(
        ' *** Error code from $SUBNAM=${INFO.i5} for M=${M.i5}, N=${N.i5}, NB=${NB.i4}, type ${IMAT.i2}');
  }

  void print9974(
      String SUBNAM, int INFO, int M, int N, int NRHS, int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM=${INFO.i5}\n ==> M =${M.i5}, N =${N.i5}, NRHS =${NRHS.i4}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9973(String SUBNAM, int INFO, int N, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5} for N =${N.i5}, type ${TYPE.i2}');
  }

  void print9972(
      String SUBNAM, int INFO, int N, int KL, int KU, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> N =${N.i5}, KL =${KL.i5}, KU =${KU.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9971(String SUBNAM, int INFO, int N, int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM=${INFO.i5} for N=${N.i5}, NB=${NB.i4}, type ${TYPE.i2}');
  }

  void print9970(String SUBNAM, int INFO, int N, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5} for N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9969(String SUBNAM, int INFO, String NORM, int N, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5} for NORM = \'${NORM.a1}\', N =${N.i5}, type ${TYPE.i2}');
  }

  void print9968(
      String SUBNAM, int INFO, String NORM, int N, int KL, int KU, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> NORM =\'${NORM.a1}\', N =${N.i5}, KL =${KL.i5}, KU =${KU.i5}, type ${TYPE.i2}');
  }

  void print9967(String SUBNAM, int INFO, String NORM, String UPLO, String DIAG,
      int N, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> NORM=\'${NORM.a1}\', UPLO =\'${UPLO.a1}\', DIAG=\'${DIAG.a1}\', N =${N.i5}, type ${TYPE.i2}');
  }

  void print9966(String SUBNAM, int INFO, String NORM, String UPLO, String DIAG,
      int N, int KD, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> NORM=\'${NORM.a1}\', UPLO =\'${UPLO.a1}\', DIAG=\'${DIAG.a1}\', N=${N.i5}, KD=${KD.i5}, type ${TYPE.i2}');
  }

  void print9965(String SUBNAM, int INFO, String TRANS, int M, int N, int NRHS,
      int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> TRANS = \'${TRANS.a1}\', M =${M.i5}, N =${N.i5}, NRHS =${NRHS.i4}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9964(String SUBNAM, int INFO, String TRANS, int N, int KL, int KU,
      int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM=${INFO.i5}\n ==> TRANS=\'${TRANS.a1}\', N =${N.i5}, KL =${KL.i5}, KU =${KU.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9963(
      String SUBNAM, int INFO, String TRANS, int N, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> TRANS = \'${TRANS.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9962(
      String SUBNAM, int INFO, String UPLO, String DIAG, int N, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO=\'${UPLO.a1}\', DIAG =\'${DIAG.a1}\', N =${N.i5}, type ${TYPE.i2}');
  }

  void print9961(String SUBNAM, int INFO, String UPLO, String DIAG, int N,
      int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO=\'${UPLO.a1}\', DIAG =\'${DIAG.a1}\', N =${N.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9960(String SUBNAM, int INFO, String UPLO, int N, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5} for UPLO = \'${UPLO.a1}\', N =${N.i5}, type ${TYPE.i2}');
  }

  void print9959(
      String SUBNAM, int INFO, String UPLO, int N, int KD, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, KD =${KD.i5}, type ${TYPE.i2}');
  }

  void print9958(
      String SUBNAM, int INFO, String UPLO, int N, int KD, int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, KD =${KD.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9957(
      String SUBNAM, int INFO, String UPLO, int N, int KD, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM=${INFO.i5}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, KD =${KD.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9956(
      String SUBNAM, int INFO, String UPLO, int N, int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }

  void print9955(
      String SUBNAM, int INFO, String UPLO, int N, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9954(String SUBNAM, int INFO, String UPLO, String TRANS,
      String DIAG, int N, int KD, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.a1}\', N=${N.i5}, KD=${KD.i5}, NRHS=${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9953(String SUBNAM, int INFO, String UPLO, String TRANS,
      String DIAG, int N, int NRHS, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.a1}\', N =${N.i5}, NRHS =${NRHS.i4}, type ${TYPE.i2}');
  }

  void print9952(String SUBNAM, int INFO, String UPLO, String TRANS,
      String DIAG, String NORMIN, int N, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.a1}\', NORMIN=\'${NORMIN.a1}\', N =${N.i5}, type ${TYPE.i2}');
  }

  void print9951(String SUBNAM, int INFO, String UPLO, String TRANS,
      String DIAG, String NORMIN, int N, int KD, int TYPE) {
    println(
        ' *** Error code from $SUBNAM =${INFO.i5}\n ==> UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.a1}\', NORMIN=\'${NORMIN.a1}\', N=${N.i5}, KD=${KD.i5}, type ${TYPE.i2}');
  }

  void print9950(String SUBNAM, int INFO) {
    println(' *** Error code from $SUBNAM =${INFO.i5}');
  }

  void print9949() {
    println(' ==> Doing only the condition estimate for this case');
  }

  void print9930(
      String SUBNAM, int INFO, int M, int N, int NX, int NB, int TYPE) {
    println(
        ' *** Error code from $SUBNAM=${INFO.i5}\n ==> M =${M.i5}, N =${N.i5}, NX =${NX.i5}, NB =${NB.i4}, type ${TYPE.i2}');
  }
}
