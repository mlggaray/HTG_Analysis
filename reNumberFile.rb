#!/usr/bin/env ruby




file = ARGV[0]
newFile = ARGV[1]

writeNewFile = File.open(newFile, 'w')



count = 1
File.open(file) { |ff|
    ff.each { |line|
      if(line =~ /^#/ )
        line.gsub!(/^#/, "")
        writeNewFile.puts line
        next
      end
    line.strip!
    writeNewFile.puts "#{count}\t#{line}"
    count += 1
  }
}

writeNewFile.close
