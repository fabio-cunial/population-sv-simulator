task tar_asm {
  File ref
  File hapOne
  File hapTwo
  String sample
  String threads
  String mem_gb
  command {
    mkdir -p asm/${sample}
    cp ${ref} asm/ref.fa
    cp ${hapOne} asm/${sample}/h1.fa.gz
    cp ${hapTwo} asm/${sample}/h2.fa.gz
    tar zcvf asm.tgz asm/
  }
  output {
    File asm_tar = "asm.tgz"
  }
}
