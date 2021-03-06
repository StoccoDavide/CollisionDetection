# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                     #
# The COLLISION DETECTION project                                     #
#                                                                     #
# Copyright (c) 2020-2021, Davide Stocco and Enrico Bertolazzi        #
#                                                                     #
# The COLLISION DETECTION project and its components are supplied     #
# under the terms of the open source BSD 2-Clause License. The        #
# contents of the COLLISION DETECTION project and its components may  #
# not be copied or disclosed except in accordance with the terms of   #
# the BSD 2-Clause License.                                           #
#                                                                     #
# URL: https://opensource.org/licenses/BSD-2-Clause                   #
#                                                                     #
#    Davide Stocco                                                    #
#    Department of Industrial Engineering                             #
#    University of Trento                                             #
#    e-mail: davide.stocco@unitn.it                                   #
#                                                                     #
#    Enrico Bertolazzi                                                #
#    Department of Industrial Engineering                             #
#    University of Trento                                             #
#    e-mail: enrico.bertolazzi@unitn.it                               #
#                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

%w(colorize fileutils pathname rubygems/package net/http zip zlib uri openssl).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require_relative "./Rakefile_conf.rb"

if (/cygwin|mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
  # LINUX
  task :default => [:install_linux]
elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
  # OSX
  task :default => [:install_osx]
else
  # WINDOWS
  task :default => [:install_windows]
end

#
# https://stackoverflow.com/questions/6934185/ruby-net-http-following-redirects/6934503
#
def url_resolve(
  uri_str,
  agent        = 'curl/7.43.0',
  max_attempts = 10,
  timeout      = 10
)
  attempts = 0
  cookie   = nil

  until attempts >= max_attempts
    attempts += 1

    url  = URI.parse(uri_str)
    http = Net::HTTP.new(url.host, url.port)
    http.open_timeout = timeout
    http.read_timeout = timeout
    path = url.path
    path = '/' if path == ''
    path += '?' + url.query unless url.query.nil?

    params = { 'User-Agent' => agent, 'Accept' => '*/*' }
    params['Cookie'] = cookie unless cookie.nil?
    request = Net::HTTP::Get.new(path, params)

    if url.instance_of?(URI::HTTPS)
      http.use_ssl = true
      http.verify_mode = OpenSSL::SSL::VERIFY_NONE
    end
    response = http.request(request)

    case response
    when Net::HTTPSuccess then
      break
    when Net::HTTPRedirection then
      location = response['Location']
      cookie   = response['Set-Cookie']
      new_uri  = URI.parse(location)
      uri_str  = if new_uri.relative? then url + location else new_uri.to_s end
    else
      raise 'Unexpected response: ' + response.inspect
    end
  end
  raise 'Too many http redirects' if attempts == max_attempts
  uri_str
  # response.body
end

def url_download( url_address, filename )
  if File.exist?(filename)
    puts "file: `#{filename}` already downloaded"
  else
    puts "downloading: #{filename}..."
    uri_str = url_resolve(url_address)
    uri     = URI( uri_str )
    File.binwrite( filename, Net::HTTP.get(uri) )
    puts "done"
  end
end

#
# https://stackoverflow.com/questions/856891/unzip-zip-tar-tag-gz-files-with-ruby
#
def extract_tgz( tar_gz_archive, destination = '.' )
  tar_longlink = '././@LongLink'
  Gem::Package::TarReader.new( Zlib::GzipReader.open tar_gz_archive ) do |tar|
    dest = nil
    tar.each do |entry|
      if entry.full_name == tar_longlink
        dest = File.join destination, entry.read.strip
        next
      end
      dest ||= File.join destination, entry.full_name
      if entry.directory?
        File.delete dest if File.file? dest
        FileUtils.mkdir_p dest, :mode => entry.header.mode, :verbose => false
      elsif entry.file?
        FileUtils.rm_rf dest if File.directory? dest
        FileUtils.mkdir_p File.dirname(dest), :mode => 0777, :verbose => false
        File.open dest, "wb" do |f| f.print entry.read end
        FileUtils.chmod entry.header.mode, dest, :verbose => false
      elsif entry.header.typeflag == '' #file?
        File.open dest, "wb" do |f| f.print entry.read end
        FileUtils.chmod entry.header.mode, dest, :verbose => false
      elsif entry.header.typeflag == '2' #Symlink!
        File.symlink entry.header.linkname, dest
      elsif entry.header.typeflag == 'g' && entry.full_name == "pax_global_header"
        puts "Skip entry: #{entry.full_name} type: #{entry.header.typeflag}."
      else
        puts "Unkown tar entry: #{entry.full_name} type: #{entry.header.typeflag}."
      end
      dest = nil
    end
  end
end

#
# https://stackoverflow.com/questions/19754883/how-to-unzip-a-zip-file-containing-folders-and-files-in-rails-while-keeping-the
#
def extract_zip( filename, destination_path='.' )
  if Zip.constants.include? :File
    zzfile = Zip::File
  else
    zzfile = Zip::ZipFile
  end
  zzfile.open(filename) do |zip_file|
    zip_file.each do |f|
      f_path=File.join(destination_path, f.name)
      FileUtils.mkdir_p(File.dirname(f_path))
      zip_file.extract(f, f_path) unless File.exist?(f_path)
    end
  end
end

def win_vs( bits, year )
  
  tmp = " -DBITS:VAR=#{bits} -DYEAR:VAR=#{year} "
  if USE_NMAKE then
    tmp = 'cmake -G "NMake Makefiles" ' + tmp
  else
    
    win32_64  = ''
    win32_64_ = '-A Win32'
    case bits
    when /x64/
      win32_64  = ' Win64'
      win32_64_ = ' -A x64'
    end

    case year
    when "2010"
      tmp = 'cmake -G "Visual Studio 10 2010' + win32_64 +'" ' + tmp
    when "2012"
      tmp = 'cmake -G "Visual Studio 11 2012' + win32_64 +'" ' + tmp
    when "2013"
      tmp = 'cmake -G "Visual Studio 12 2013' + win32_64 +'" ' + tmp
    when "2015"
      tmp = 'cmake -G "Visual Studio 14 2015' + win32_64 +'" ' + tmp
    when "2017"
      tmp = 'cmake -G "Visual Studio 15 2017' + win32_64 +'" ' + tmp
    when "2019"
      tmp = 'cmake -G "Visual Studio 16 2019"' + win32_64_ + tmp
    else
      puts "Visual Studio year #{year} not supported!\n";
      return ""
    end
  end
  return tmp
end
