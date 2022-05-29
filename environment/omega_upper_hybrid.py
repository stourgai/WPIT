def omega_upper_hybrid(wce_arg,wpe_arg,wci_arg,wpi_arg):
    fac1=(wce_arg*wce_arg+wpe_arg*wpe_arg+wci_arg*wci_arg+wpi_arg*wpi_arg)/2
    fac2=(wce_arg*wce_arg+wpe_arg*wpe_arg-wci_arg*wci_arg-wpi_arg*wpi_arg)**2+4*wpe_arg*wpe_arg*wpi_arg*wpi_arg
    wuh_tmp=fac1+0.5*fac2
    return wuh_tmp
