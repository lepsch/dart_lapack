      SUBROUTINE SLAHILB( N, NRHS, A, LDA, X, LDX, B, LDB, WORK, INFO);

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     N, NRHS, LDA, LDX, LDB, INFO;
      // .. Array Arguments ..
      REAL A(LDA, N), X(LDX, NRHS), B(LDB, NRHS), WORK(N);
      // ..

*  =====================================================================
      // .. Local Scalars ..
      int     TM, TI, R;
      int     M;
      int     I, J;
      // ..
      // .. Parameters ..
      // NMAX_EXACT   the largest dimension where the generated data is
                   // exact.
      // NMAX_APPROX  the largest dimension where the generated data has
                   // a small componentwise relative error.
      int     NMAX_EXACT, NMAX_APPROX;
      const     NMAX_EXACT = 6, NMAX_APPROX = 11;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. External Functions
      // EXTERNAL SLASET
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      if (N < 0 || N > NMAX_APPROX) {
         INFO = -1;
      } else if (NRHS < 0) {
         INFO = -2;
      } else if (LDA < N) {
         INFO = -4;
      } else if (LDX < N) {
         INFO = -6;
      } else if (LDB < N) {
         INFO = -8;
      }
      if (INFO < 0) {
         xerbla('SLAHILB', -INFO);
         RETURN;
      }
      if (N > NMAX_EXACT) {
         INFO = 1;
      }

      // Compute M = the LCM of the integers [1, 2*N-1].  The largest
      // reasonable N is small enough that integers suffice (up to N = 11).
      M = 1;
      for (I = 2; I <= (2*N-1); I++) {
         TM = M;
         TI = I;
         R = MOD(TM, TI);
         DO WHILE (R != 0);
            TM = TI;
            TI = R;
            R = MOD(TM, TI);
         }
         M = (M / TI) * I;
      }

      // Generate the scaled Hilbert matrix in A
      for (J = 1; J <= N; J++) {
         for (I = 1; I <= N; I++) {
            A(I, J) = REAL(M) / (I + J - 1);
         }
      }

      // Generate matrix B as simply the first NRHS columns of M * the
      // identity.
      slaset('Full', N, NRHS, 0.0, REAL(M), B, LDB);

      // Generate the true solutions in X.  Because B = the first NRHS
      // columns of M*I, the true solutions are just the first NRHS columns
      // of the inverse Hilbert matrix.
      WORK(1) = N;
      for (J = 2; J <= N; J++) {
         WORK(J) = (  ( (WORK(J-1)/(J-1)) * (J-1 - N) ) /(J-1)  ) * (N +J -1);
      }

      for (J = 1; J <= NRHS; J++) {
         for (I = 1; I <= N; I++) {
            X(I, J) = (WORK(I)*WORK(J)) / (I + J - 1);
         }
      }

      }
