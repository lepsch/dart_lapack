void slags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  bool UPPER;
  double A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, SNQ, SNU, SNV;
  // ..

// =====================================================================

  // .. Parameters ..
  double ZERO;
  const ZERO = 0.0;
  // ..
  // .. Local Scalars ..
  double A,
      AUA11,
      AUA12,
      AUA21,
      AUA22,
      AVB11,
      AVB12,
      AVB21,
      AVB22,
      CSL,
      CSR,
      D,
      S1,
      S2,
      SNL,
      SNR,
      UA11R,
      UA22R,
      VB11R,
      VB22R,
      B,
      C,
      R,
      UA11,
      UA12,
      UA21,
      UA22,
      VB11,
      VB12,
      VB21,
      VB22;
  // ..
  // .. External Subroutines ..
  // EXTERNAL SLARTG, SLASV2
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC ABS
  // ..
  // .. Executable Statements ..

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

    slasv2(A, B, D, S1, S2, SNR, CSR, SNL, CSL);

    if ((CSL).abs() >= (SNL).abs() || (CSR).abs() >= (SNR).abs()) {
      // Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
      // and (1,2) element of |U|**T *|A| and |V|**T *|B|.

      UA11R = CSL * A1;
      UA12 = CSL * A2 + SNL * A3;

      VB11R = CSR * B1;
      VB12 = CSR * B2 + SNR * B3;

      AUA12 = (CSL).abs() * (A2).abs() + (SNL).abs() * (A3).abs();
      AVB12 = (CSR).abs() * (B2).abs() + (SNR).abs() * (B3).abs();

      // zero (1,2) elements of U**T *A and V**T *B

      if (((UA11R).abs() + (UA12).abs()) != ZERO) {
        if (AUA12 / ((UA11R).abs() + (UA12).abs()) <=
            AVB12 / ((VB11R).abs() + (VB12).abs())) {
          slartg(-UA11R, UA12, CSQ, SNQ, R);
        } else {
          slartg(-VB11R, VB12, CSQ, SNQ, R);
        }
      } else {
        slartg(-VB11R, VB12, CSQ, SNQ, R);
      }

      CSU = CSL;
      SNU = -SNL;
      CSV = CSR;
      SNV = -SNR;
    } else {
      // Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
      // and (2,2) element of |U|**T *|A| and |V|**T *|B|.

      UA21 = -SNL * A1;
      UA22 = -SNL * A2 + CSL * A3;

      VB21 = -SNR * B1;
      VB22 = -SNR * B2 + CSR * B3;

      AUA22 = (SNL).abs() * (A2).abs() + (CSL).abs() * (A3).abs();
      AVB22 = (SNR).abs() * (B2).abs() + (CSR).abs() * (B3).abs();

      // zero (2,2) elements of U**T*A and V**T*B, and then swap.

      if (((UA21).abs() + (UA22).abs()) != ZERO) {
        if (AUA22 / ((UA21).abs() + (UA22).abs()) <=
            AVB22 / ((VB21).abs() + (VB22).abs())) {
          slartg(-UA21, UA22, CSQ, SNQ, R);
        } else {
          slartg(-VB21, VB22, CSQ, SNQ, R);
        }
      } else {
        slartg(-VB21, VB22, CSQ, SNQ, R);
      }

      CSU = SNL;
      SNU = CSL;
      CSV = SNR;
      SNV = CSR;
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

    slasv2(A, C, D, S1, S2, SNR, CSR, SNL, CSL);

    if ((CSR).abs() >= (SNR).abs() || (CSL).abs() >= (SNL).abs()) {
      // Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
      // and (2,1) element of |U|**T *|A| and |V|**T *|B|.

      UA21 = -SNR * A1 + CSR * A2;
      UA22R = CSR * A3;

      VB21 = -SNL * B1 + CSL * B2;
      VB22R = CSL * B3;

      AUA21 = (SNR).abs() * (A1).abs() + (CSR).abs() * (A2).abs();
      AVB21 = (SNL).abs() * (B1).abs() + (CSL).abs() * (B2).abs();

      // zero (2,1) elements of U**T *A and V**T *B.

      if (((UA21).abs() + (UA22R).abs()) != ZERO) {
        if (AUA21 / ((UA21).abs() + (UA22R).abs()) <=
            AVB21 / ((VB21).abs() + (VB22R).abs())) {
          slartg(UA22R, UA21, CSQ, SNQ, R);
        } else {
          slartg(VB22R, VB21, CSQ, SNQ, R);
        }
      } else {
        slartg(VB22R, VB21, CSQ, SNQ, R);
      }

      CSU = CSR;
      SNU = -SNR;
      CSV = CSL;
      SNV = -SNL;
    } else {
      // Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
      // and (1,1) element of |U|**T *|A| and |V|**T *|B|.

      UA11 = CSR * A1 + SNR * A2;
      UA12 = SNR * A3;

      VB11 = CSL * B1 + SNL * B2;
      VB12 = SNL * B3;

      AUA11 = (CSR).abs() * (A1).abs() + (SNR).abs() * (A2).abs();
      AVB11 = (CSL).abs() * (B1).abs() + (SNL).abs() * (B2).abs();

      // zero (1,1) elements of U**T*A and V**T*B, and then swap.

      if (((UA11).abs() + (UA12).abs()) != ZERO) {
        if (AUA11 / ((UA11).abs() + (UA12).abs()) <=
            AVB11 / ((VB11).abs() + (VB12).abs())) {
          slartg(UA12, UA11, CSQ, SNQ, R);
        } else {
          slartg(VB12, VB11, CSQ, SNQ, R);
        }
      } else {
        slartg(VB12, VB11, CSQ, SNQ, R);
      }

      CSU = SNR;
      SNU = CSR;
      CSV = SNL;
      SNV = CSL;
    }
  }

  return;
}
