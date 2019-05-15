procs = [1,2,4,8]
iter = []
procs.each do |proc|
  puts "============Start Processors: #{proc}==================="
  accuracy_criteria = 0.1
  while accuracy_criteria >= 0.000001
    start_time = Time.now
    op = `mpirun -np #{proc} ./monte_carlo_pi #{accuracy_criteria}`
    end_time = Time.now
    iter << [accuracy_criteria, end_time - start_time, proc]
    accuracy_criteria = accuracy_criteria/10
  end
  puts "========================================================"
end

require 'csv'
CSV.open("monte_carlo_pi.csv", "w") do |csv|
  csv << ['Accuracy Criteria', 'Time in Seconds', 'Processors']
  iter.each do |x|
    csv << x
  end
end
