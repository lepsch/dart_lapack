import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlasv2.dart';

void dlags2(
  final bool UPPER,
  final double A1,
  final double A2,
  final double A3,
  final double B1,
  final double B2,
  final double B3,
  final Box<double> CSU,
  final Box<double> SNU,
  final Box<double> CSV,
  final Box<double> SNV,
  final Box<double> CSQ,
  final Box<double> SNQ,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  double A,
      AUA11,
      AUA12,
      AUA21,
      AUA22,
      AVB11,
      AVB12,
      AVB21,
      AVB22,
      B,
      C,
      D,
      UA11,
      UA11R,
      UA12,
      UA21,
      UA22,
      UA22R,
      VB11,
      VB11R,
      VB12,
      VB21,
      VB22,
      VB22R;
  final CSL = Box(0.0),
      CSR = Box(0.0),
      S1 = Box(0.0),
      S2 = Box(0.0),
      R = Box(0.0),
      SNL = Box(0.0),
      SNR = Box(0.0);

  if (UPPER) {
    // Input matrices A and B are upper triangular matrices

    // Form matrix C = A*adj(B) = ( a b )
    // ( 0 d )

    A = A1 * B3;
    D = A3 * B1;
    B = A2 * B1 - A1 * B2;

    // The SVD of real 2-by-2 triangular C

    // ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
    // ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )

    dlasv2(A, B, D, S1, S2, SNR, CSR, SNL, CSL);

    if (CSL.value.abs() >= SNL.value.abs() ||
        CSR.value.abs() >= SNR.value.abs()) {
      // Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
      // and (1,2) element of |U|**T *|A| and |V|**T *|B|.

      UA11R = CSL.value * A1;
      UA12 = CSL.value * A2 + SNL.value * A3;

      VB11R = CSR.value * B1;
      VB12 = CSR.value * B2 + SNR.value * B3;

      AUA12 = CSL.value.abs() * A2.abs() + SNL.value.abs() * A3.abs();
      AVB12 = CSR.value.abs() * B2.abs() + SNR.value.abs() * B3.abs();

      // zero (1,2) elements of U**T *A and V**T *B

      if ((UA11R.abs() + UA12.abs()) != ZERO) {
        if (AUA12 / (UA11R.abs() + UA12.abs()) <=
            AVB12 / (VB11R.abs() + VB12.abs())) {
          dlartg(-UA11R, UA12, CSQ, SNQ, R);
        } else {
          dlartg(-VB11R, VB12, CSQ, SNQ, R);
        }
      } else {
        dlartg(-VB11R, VB12, CSQ, SNQ, R);
      }

      CSU.value = CSL.value;
      SNU.value = -SNL.value;
      CSV.value = CSR.value;
      SNV.value = -SNR.value;
    } else {
      // Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
      // and (2,2) element of |U|**T *|A| and |V|**T *|B|.

      UA21 = -SNL.value * A1;
      UA22 = -SNL.value * A2 + CSL.value * A3;

      VB21 = -SNR.value * B1;
      VB22 = -SNR.value * B2 + CSR.value * B3;

      AUA22 = SNL.value.abs() * A2.abs() + CSL.value.abs() * A3.abs();
      AVB22 = SNR.value.abs() * B2.abs() + CSR.value.abs() * B3.abs();

      // zero (2,2) elements of U**T*A and V**T*B, and then swap.

      if ((UA21.abs() + UA22.abs()) != ZERO) {
        if (AUA22 / (UA21.abs() + UA22.abs()) <=
            AVB22 / (VB21.abs() + VB22.abs())) {
          dlartg(-UA21, UA22, CSQ, SNQ, R);
        } else {
          dlartg(-VB21, VB22, CSQ, SNQ, R);
        }
      } else {
        dlartg(-VB21, VB22, CSQ, SNQ, R);
      }

      CSU.value = SNL.value;
      SNU.value = CSL.value;
      CSV.value = SNR.value;
      SNV.value = CSR.value;
    }
  } else {
    // Input matrices A and B are lower triangular matrices

    // Form matrix C = A*adj(B) = ( a 0 )
    // ( c d )

    A = A1 * B3;
    D = A3 * B1;
    C = A2 * B3 - A3 * B2;

    // The SVD of real 2-by-2 triangular C

    // ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
    // ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )

    dlasv2(A, C, D, S1, S2, SNR, CSR, SNL, CSL);

    if (CSR.value.abs() >= SNR.value.abs() ||
        CSL.value.abs() >= SNL.value.abs()) {
      // Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
      // and (2,1) element of |U|**T *|A| and |V|**T *|B|.

      UA21 = -SNR.value * A1 + CSR.value * A2;
      UA22R = CSR.value * A3;

      VB21 = -SNL.value * B1 + CSL.value * B2;
      VB22R = CSL.value * B3;

      AUA21 = SNR.value.abs() * A1.abs() + CSR.value.abs() * A2.abs();
      AVB21 = SNL.value.abs() * B1.abs() + CSL.value.abs() * B2.abs();

      // zero (2,1) elements of U**T *A and V**T *B.

      if ((UA21.abs() + UA22R.abs()) != ZERO) {
        if (AUA21 / (UA21.abs() + UA22R.abs()) <=
            AVB21 / (VB21.abs() + VB22R.abs())) {
          dlartg(UA22R, UA21, CSQ, SNQ, R);
        } else {
          dlartg(VB22R, VB21, CSQ, SNQ, R);
        }
      } else {
        dlartg(VB22R, VB21, CSQ, SNQ, R);
      }

      CSU.value = CSR.value;
      SNU.value = -SNR.value;
      CSV.value = CSL.value;
      SNV.value = -SNL.value;
    } else {
      // Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
      // and (1,1) element of |U|**T *|A| and |V|**T *|B|.

      UA11 = CSR.value * A1 + SNR.value * A2;
      UA12 = SNR.value * A3;

      VB11 = CSL.value * B1 + SNL.value * B2;
      VB12 = SNL.value * B3;

      AUA11 = CSR.value.abs() * A1.abs() + SNR.value.abs() * A2.abs();
      AVB11 = CSL.value.abs() * B1.abs() + SNL.value.abs() * B2.abs();

      // zero (1,1) elements of U**T*A and V**T*B, and then swap.

      if ((UA11.abs() + UA12.abs()) != ZERO) {
        if (AUA11 / (UA11.abs() + UA12.abs()) <=
            AVB11 / (VB11.abs() + VB12.abs())) {
          dlartg(UA12, UA11, CSQ, SNQ, R);
        } else {
          dlartg(VB12, VB11, CSQ, SNQ, R);
        }
      } else {
        dlartg(VB12, VB11, CSQ, SNQ, R);
      }

      CSU.value = SNR.value;
      SNU.value = CSR.value;
      CSV.value = SNL.value;
      SNV.value = CSL.value;
    }
  }
}
