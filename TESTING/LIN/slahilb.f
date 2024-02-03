      SUBROUTINE SLAHILB( N, NRHS, A, LDA, X, LDX, B, LDB, WORK, INFO)

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     N, NRHS, LDA, LDX, LDB, INFO;
      // .. Array Arguments ..
      REAL A(LDA, N), X(LDX, NRHS), B(LDB, NRHS), WORK(N)
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
      // .. External Functions
      // EXTERNAL SLASET
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if (N .LT. 0 .OR. N .GT. NMAX_APPROX) {
         INFO = -1
      } else if (NRHS .LT. 0) {
         INFO = -2
      } else if (LDA .LT. N) {
         INFO = -4
      } else if (LDX .LT. N) {
         INFO = -6
      } else if (LDB .LT. N) {
         INFO = -8
      }
      if (INFO .LT. 0) {
         CALL XERBLA('SLAHILB', -INFO)
         RETURN
      }
      if (N .GT. NMAX_EXACT) {
         INFO = 1
      }

      // Compute M = the LCM of the integers [1, 2*N-1].  The largest
      // reasonable N is small enough that integers suffice (up to N = 11).
      M = 1
      DO I = 2, (2*N-1)
         TM = M
         TI = I
         R = MOD(TM, TI)
         DO WHILE (R .NE. 0)
            TM = TI
            TI = R
            R = MOD(TM, TI)
         END DO
         M = (M / TI) * I
      END DO

      // Generate the scaled Hilbert matrix in A
      DO J = 1, N
         DO I = 1, N
            A(I, J) = REAL(M) / (I + J - 1)
         END DO
      END DO

      // Generate matrix B as simply the first NRHS columns of M * the
      // identity.
      CALL SLASET('Full', N, NRHS, 0.0, REAL(M), B, LDB)

      // Generate the true solutions in X.  Because B = the first NRHS
      // columns of M*I, the true solutions are just the first NRHS columns
      // of the inverse Hilbert matrix.
      WORK(1) = N
      DO J = 2, N
         WORK(J) = (  ( (WORK(J-1)/(J-1)) * (J-1 - N) ) /(J-1)  ) * (N +J -1)
      END DO

      DO J = 1, NRHS
         DO I = 1, N
            X(I, J) = (WORK(I)*WORK(J)) / (I + J - 1)
         END DO
      END DO

      }
