# Introduction

This code will generate distributed multi-channel simulated speech data based on librispeech corpus.

The simulation environment is that a single speaker is located in a reverberant environment, containing point noise and background noise, and the microphones are randomly placed in the room.

You can modify the number of channels, room size, reverberation time(T60), noise type and signal-to-noise ratio to obtain the desired simulation data. You can

The File structure generation mechanism refers to [DIRHA_English_phrich](https://github.com/SHINE-FBK/DIRHA_English_phrich)

Use [RIR-Generator](https://github.com/ehabets/RIR-Generator) to generate room impulse response

# Generate default simulated data

0. Download Librispeech data set form here: [LibriSpeech ASR corpus](http://www.openslr.org/12/)

1. Open Matlab and Open "**generate/main.m**"

2. Set "**debug_mode**" to `false`

3. Set "**librispeech_dir**" to the path where you store data, 

   for example: `librispeech_dir = '/home/public/data/LibriSpeech'`;

4. Set "**noise_dir**", "**diffuse_noise_dir**" and "**point_noise_dir**" to the path where you store noise data, For the specific settings of these folders, see **Generate specific simulated data**. You need to prepare the noise data yourself.

5. Set "**setname**" to the sets you expect to generate, 

   for example: `setname = {'train-clean-100'}; ` or `setname = {'test-clean', 'test-other'}; `

6. Set "**target_dir**" to the path you expect to output , 

   for example: `target_dir = '/home/public/Generate_data/LibriSpeech';`

7. Run "**generate/main.m**"



# Output data structure

After the data is generated, the output folder structure is as follows:

> target_dir
>
> >train-clean-100
> >
> >> noise_speech
> >>
> >> direct_speech
> >>
> >> bkg_noise
> >>
> >> info
> >
> >other set
> >
> >>...

The file structure of "**noise_speech**", "**direct_speech**" and "**bkg_noise**" directory is the same as the original data. 

- **noise_speech**: contains noisy speech data
- **direct_speech**: contains direct speech data without noise and reverberation 
- **bkg_noise**: contains background noise data

**noise_speech data = direct_speech data + bkg_noise data**

You can view the information of the generated data in "**info**" directory



# Generate specific simulated data

* In "**generate/main.m**", you can modify the following parameters:

  * "**mic**"

    Default is `{'ad-hoc',16}` , means ad-hoc array with 16 channels, you can change the number of channels by modifying this parameter to `{'ad-hoc',30}`

  * "**option.value**"

    Default is `[1,1,1,1]`, corresponds to "option.key", set to 1 means True, set to 0 means False

  * "**max_pum**"

    Default is `4`, means the numbers of point source noise is randomly selected from 1 to 4

  * "**noise_dir**", "**diffuse_noise_dir**" and "**point_noise_dir**"

    The diffuse noise and point noise in the training set come from "noise_dir"

    The diffuse noise in the test set comes from "diffuse_noise_dir"

    The point noise in the test set comes from "point_noise_dir"

    

* In "**generate/noise_speech_generate.m**"

  * "**snr_range_for_diffuse_noise**"

    SNR range of original speech and original diffuse noise

  * "**snr_range_for_point_noise**"

    SNR range of original speech and original point noise

  * "**room_shape**"

    The range of room shape

  * "**min_distance_source2wall**"

    The minimum distance between source to wall

  * "**min_distance_source2mic**"

    The minimum distance between source to microphone

  * "**T60**", "**T60_sample_mean**" and "**T60_sample_std**"

    The range of T60. The mean and  standard deviation of T60 sampling.
    
    

# Generating the WSJ0-adhoc dataset for speech separation

1. You can get n-speakers-mixed dataset for speech separation through regarding inference speakers as point noise.

2. Set "**librispeech_dir**" to the path where you store row data which can be librispeech or WSJ0.

3. You can get the file names of different speakers to generate mixed speech by files under folder `.../file_name/`. It concludes three files which are used to  generate train (tr), validation (cv) and test (tt) datasets respectively. 

   In each file, you can get four information: the file name of speaker 1, the file name of speaker 2, the SNR of speaker 1 and speaker 2 and the file name of auxiliary speaker.

4. You can set the remaining parameters according to your needs.




​    
