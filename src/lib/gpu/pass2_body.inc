! One cell of the fused flux pass 2 (Flux_Pass2_Fused): the six face
! buffers summed into R_flat(c0,:) in a fixed order; shared so the box and the list sweeps match.
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nprim; tmp(s) = 0d0; end do
      ! convective i (+ right face, - left face)
      if (i < n1) then
        f = (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
        do s = 1, nprim; tmp(s) = tmp(s) + fhat(f, s); end do
      end if
      if (i > 1) then
        f = (i-2) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
        do s = 1, nprim; tmp(s) = tmp(s) - fhat(f, s); end do
      end if
      ! convective j
      if (j < n2) then
        f = (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
        do s = 1, nprim; tmp(s) = tmp(s) + fhat_j(f, s); end do
      end if
      if (j > 1) then
        f = (i-1) + (j-2)*n1 + (k-1)*n1*(n2-1) + 1
        do s = 1, nprim; tmp(s) = tmp(s) - fhat_j(f, s); end do
      end if
      ! convective k
      if (k < n3) then
        f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
        do s = 1, nprim; tmp(s) = tmp(s) + fhat_k(f, s); end do
      end if
      if (k > 1) then
        f = (i-1) + (j-1)*n1 + (k-2)*n1*n2 + 1
        do s = 1, nprim; tmp(s) = tmp(s) - fhat_k(f, s); end do
      end if
      if (visc_l) then
        ! diffusive i (- right face, + left face)
        if (i < n1) then
          f = (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
          do s = 1, nprim; tmp(s) = tmp(s) - dfhat(f, s); end do
        end if
        if (i > 1) then
          f = (i-2) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
          do s = 1, nprim; tmp(s) = tmp(s) + dfhat(f, s); end do
        end if
        ! diffusive j
        if (j < n2) then
          f = (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
          do s = 1, nprim; tmp(s) = tmp(s) - dfhat_j(f, s); end do
        end if
        if (j > 1) then
          f = (i-1) + (j-2)*n1 + (k-1)*n1*(n2-1) + 1
          do s = 1, nprim; tmp(s) = tmp(s) + dfhat_j(f, s); end do
        end if
        ! diffusive k
        if (k < n3) then
          f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
          do s = 1, nprim; tmp(s) = tmp(s) - dfhat_k(f, s); end do
        end if
        if (k > 1) then
          f = (i-1) + (j-1)*n1 + (k-2)*n1*n2 + 1
          do s = 1, nprim; tmp(s) = tmp(s) + dfhat_k(f, s); end do
        end if
      end if
      do s = 1, nprim; R_flat(c0,s) = tmp(s); end do
