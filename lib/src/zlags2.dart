import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlasv2.dart';
import 'package:lapack/src/zlartg.dart';

void zlags2(
  final bool UPPER,
  final double A1,
  final Complex A2,
  final double A3,
  final double B1,
  final Complex B2,
  final double B3,
  final Box<double> CSU,
  final Box<Complex> SNU,
  final Box<double> CSV,
  final Box<Complex> SNV,
  final Box<double> CSQ,
  final Box<Complex> SNQ,
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
      AVB12,
      AVB11,
      AVB21,
      AVB22,
      D,
      FB,
      FC,
      UA11R,
      UA22R,
      VB11R,
      VB22R;
  Complex B, C, D1, UA11, UA12, UA21, UA22, VB11, VB12, VB21, VB22;
  final S1 = Box(0.0),
      S2 = Box(0.0),
      SNR = Box(0.0),
      CSR = Box(0.0),
      SNL = Box(0.0),
      CSL = Box(0.0);
  final R = Box(Complex.zero);

  double ABS1(Complex T) => T.toDouble().abs() + T.imaginary.abs();

  if (UPPER) {
    // Input matrices A and B are upper triangular matrices

    // Form matrix C = A*adj(B) = ( a b )
    //                            ( 0 d )

    A = A1 * B3;
    D = A3 * B1;
    B = A2 * B1.toComplex() - A1.toComplex() * B2;
    FB = (B).abs();

    // Transform complex 2-by-2 matrix C to real matrix by unitary
    // diagonal matrix diag(1,D1).

    D1 = Complex.one;
    if (FB != ZERO) D1 = B / FB.toComplex();

    // The SVD of real 2 by 2 triangular C

    // ( CSL.value -SNL.value )*( A B )*(  CSR.value  SNR.value ) = ( R.value 0 )
    // ( SNL.value  CSL.value ) ( 0 D ) ( -SNR.value  CSR.value )   ( 0 T )

    dlasv2(A, FB, D, S1, S2, SNR, CSR, SNL, CSL);

    if ((CSL.value).abs() >= (SNL.value).abs() ||
        (CSR.value).abs() >= (SNR.value).abs()) {
      // Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
      // and (1,2) element of |U|**H *|A| and |V|**H *|B|.

      UA11R = CSL.value * A1;
      UA12 = CSL.value.toComplex() * A2 +
          D1 * SNL.value.toComplex() * A3.toComplex();

      VB11R = CSR.value * B1;
      VB12 = CSR.value.toComplex() * B2 +
          D1 * SNR.value.toComplex() * B3.toComplex();

      AUA12 = (CSL.value).abs() * ABS1(A2) + (SNL.value).abs() * (A3).abs();
      AVB12 = (CSR.value).abs() * ABS1(B2) + (SNR.value).abs() * (B3).abs();

      // zero (1,2) elements of U**H *A and V**H *B

      if (((UA11R).abs() + ABS1(UA12)) == ZERO) {
        zlartg(-VB11R.toComplex(), VB12.conjugate(), CSQ, SNQ, R);
      } else if (((VB11R).abs() + ABS1(VB12)) == ZERO) {
        zlartg(-UA11R.toComplex(), UA12.conjugate(), CSQ, SNQ, R);
      } else if (AUA12 / ((UA11R).abs() + ABS1(UA12)) <=
          AVB12 / ((VB11R).abs() + ABS1(VB12))) {
        zlartg(-UA11R.toComplex(), UA12.conjugate(), CSQ, SNQ, R);
      } else {
        zlartg(-VB11R.toComplex(), VB12.conjugate(), CSQ, SNQ, R);
      }

      CSU.value = CSL.value;
      SNU.value = -D1 * SNL.value.toComplex();
      CSV.value = CSR.value;
      SNV.value = -D1 * SNR.value.toComplex();
    } else {
      // Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
      // and (2,2) element of |U|**H *|A| and |V|**H *|B|.

      UA21 = -D1.conjugate() * (SNL.value * A1).toComplex();
      UA22 = -D1.conjugate() * SNL.value.toComplex() * A2 +
          (CSL.value * A3).toComplex();

      VB21 = -D1.conjugate() * (SNR.value * B1).toComplex();
      VB22 = -D1.conjugate() * SNR.value.toComplex() * B2 +
          (CSR.value * B3).toComplex();

      AUA22 = (SNL.value).abs() * ABS1(A2) + (CSL.value).abs() * (A3).abs();
      AVB22 = (SNR.value).abs() * ABS1(B2) + (CSR.value).abs() * (B3).abs();

      // zero (2,2) elements of U**H *A and V**H *B, and then swap.

      if ((ABS1(UA21) + ABS1(UA22)) == ZERO) {
        zlartg(-VB21.conjugate(), VB22.conjugate(), CSQ, SNQ, R);
      } else if ((ABS1(VB21) + (VB22).abs()) == ZERO) {
        zlartg(-UA21.conjugate(), UA22.conjugate(), CSQ, SNQ, R);
      } else if (AUA22 / (ABS1(UA21) + ABS1(UA22)) <=
          AVB22 / (ABS1(VB21) + ABS1(VB22))) {
        zlartg(-UA21.conjugate(), UA22.conjugate(), CSQ, SNQ, R);
      } else {
        zlartg(-VB21.conjugate(), VB22.conjugate(), CSQ, SNQ, R);
      }

      CSU.value = SNL.value;
      SNU.value = D1 * CSL.value.toComplex();
      CSV.value = SNR.value;
      SNV.value = D1 * CSR.value.toComplex();
    }
  } else {
    // Input matrices A and B are lower triangular matrices

    // Form matrix C = A*adj(B) = ( a 0 )
    //                            ( c d )

    A = A1 * B3;
    D = A3 * B1;
    C = A2 * B3.toComplex() - A3.toComplex() * B2;
    FC = (C).abs();

    // Transform complex 2-by-2 matrix C to real matrix by unitary
    // diagonal matrix diag(d1,1).

    D1 = Complex.one;
    if (FC != ZERO) D1 = C / FC.toComplex();

    // The SVD of real 2 by 2 triangular C

    // ( CSL.value -SNL.value )*( A 0 )*(  CSR.value  SNR.value ) = ( R.value 0 )
    // ( SNL.value  CSL.value ) ( C D ) ( -SNR.value  CSR.value )   ( 0 T )

    dlasv2(A, FC, D, S1, S2, SNR, CSR, SNL, CSL);

    if ((CSR.value).abs() >= (SNR.value).abs() ||
        (CSL.value).abs() >= (SNL.value).abs()) {
      // Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
      // and (2,1) element of |U|**H *|A| and |V|**H *|B|.

      UA21 = -D1 * (SNR.value * A1).toComplex() + CSR.value.toComplex() * A2;
      UA22R = CSR.value * A3;

      VB21 = -D1 * (SNL.value * B1).toComplex() + CSL.value.toComplex() * B2;
      VB22R = CSL.value * B3;

      AUA21 = (SNR.value).abs() * (A1).abs() + (CSR.value).abs() * ABS1(A2);
      AVB21 = (SNL.value).abs() * (B1).abs() + (CSL.value).abs() * ABS1(B2);

      // zero (2,1) elements of U**H *A and V**H *B.

      if ((ABS1(UA21) + (UA22R).abs()) == ZERO) {
        zlartg(VB22R.toComplex(), VB21, CSQ, SNQ, R);
      } else if ((ABS1(VB21) + (VB22R).abs()) == ZERO) {
        zlartg(UA22R.toComplex(), UA21, CSQ, SNQ, R);
      } else if (AUA21 / (ABS1(UA21) + (UA22R).abs()) <=
          AVB21 / (ABS1(VB21) + (VB22R).abs())) {
        zlartg(UA22R.toComplex(), UA21, CSQ, SNQ, R);
      } else {
        zlartg(VB22R.toComplex(), VB21, CSQ, SNQ, R);
      }

      CSU.value = CSR.value;
      SNU.value = -D1.conjugate() * SNR.value.toComplex();
      CSV.value = CSL.value;
      SNV.value = -D1.conjugate() * SNL.value.toComplex();
    } else {
      // Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
      // and (1,1) element of |U|**H *|A| and |V|**H *|B|.

      UA11 = (CSR.value * A1).toComplex() +
          D1.conjugate() * SNR.value.toComplex() * A2;
      UA12 = D1.conjugate() * (SNR.value * A3).toComplex();

      VB11 = (CSL.value * B1).toComplex() +
          D1.conjugate() * SNL.value.toComplex() * B2;
      VB12 = D1.conjugate() * (SNL.value * B3).toComplex();

      AUA11 = (CSR.value).abs() * (A1).abs() + (SNR.value).abs() * ABS1(A2);
      AVB11 = (CSL.value).abs() * (B1).abs() + (SNL.value).abs() * ABS1(B2);

      // zero (1,1) elements of U**H *A and V**H *B, and then swap.

      if ((ABS1(UA11) + ABS1(UA12)) == ZERO) {
        zlartg(VB12, VB11, CSQ, SNQ, R);
      } else if ((ABS1(VB11) + ABS1(VB12)) == ZERO) {
        zlartg(UA12, UA11, CSQ, SNQ, R);
      } else if (AUA11 / (ABS1(UA11) + ABS1(UA12)) <=
          AVB11 / (ABS1(VB11) + ABS1(VB12))) {
        zlartg(UA12, UA11, CSQ, SNQ, R);
      } else {
        zlartg(VB12, VB11, CSQ, SNQ, R);
      }

      CSU.value = SNR.value;
      SNU.value = D1.conjugate() * CSR.value.toComplex();
      CSV.value = SNL.value;
      SNV.value = D1.conjugate() * CSL.value.toComplex();
    }
  }
}
