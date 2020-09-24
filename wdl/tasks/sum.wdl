version 1.0

workflow runSum {
    call sum
}

task sum {
    input {
        Array[Int] integers
    }

    command <<<
        echo $((~{sep="+" integers}))
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        docker: "tpesout/hpp_base:latest"
    }
}