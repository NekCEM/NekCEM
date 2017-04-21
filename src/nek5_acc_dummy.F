 
      logical function acc_nek_present(a,n)
#ifdef _OPENACC
      use openacc
#endif
      integer n
      real    a(n)
#ifdef _OPENACC
      acc_nek_present = acc_is_present(a,n)
#else
      acc_nek_present = .false.
#endif
      return
      end
