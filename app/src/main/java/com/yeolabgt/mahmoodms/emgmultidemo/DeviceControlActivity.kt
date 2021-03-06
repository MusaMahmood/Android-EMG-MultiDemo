package com.yeolabgt.mahmoodms.emgmultidemo

import android.annotation.SuppressLint
import android.app.Activity
import android.app.ProgressDialog
import android.bluetooth.*
import android.content.Context
import android.content.DialogInterface
import android.content.Intent
import android.content.pm.ActivityInfo
import android.graphics.Color
import android.graphics.Typeface
import android.graphics.drawable.ColorDrawable
import android.media.MediaPlayer
import android.net.Uri
import android.os.Build
import android.os.Bundle
import android.os.Handler
import android.support.v4.app.NavUtils
import android.support.v4.content.FileProvider
import android.util.Log
import android.view.*
import android.widget.Button
import android.widget.TextView
import android.widget.Toast
import android.widget.ToggleButton
import com.androidplot.util.Redrawer
import com.google.common.primitives.Doubles
import com.parrot.arsdk.arcommands.ARCOMMANDS_JUMPINGSUMO_MEDIARECORDEVENT_PICTUREEVENTCHANGED_ERROR_ENUM
import com.parrot.arsdk.arcommands.ARCOMMANDS_MINIDRONE_MEDIARECORDEVENT_PICTUREEVENTCHANGED_ERROR_ENUM
import com.parrot.arsdk.arcommands.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM
import com.parrot.arsdk.arcontroller.ARCONTROLLER_DEVICE_STATE_ENUM
import com.parrot.arsdk.arcontroller.ARControllerCodec
import com.parrot.arsdk.arcontroller.ARFrame
import com.parrot.arsdk.ardiscovery.ARDiscoveryDeviceService
import com.yeolabgt.mahmoodms.actblelibrary.ActBle
import com.yeolabgt.mahmoodms.emgmultidemo.ParrotDrone.JSDrone
import com.yeolabgt.mahmoodms.emgmultidemo.ParrotDrone.MiniDrone
import kotlinx.android.synthetic.main.activit_dev_ctrl_alt.*
import org.tensorflow.contrib.android.TensorFlowInferenceInterface
import java.io.IOException
import java.text.SimpleDateFormat
import java.util.*
import kotlin.experimental.and
import kotlin.experimental.or

/**
 * Created by mahmoodms on 5/31/2016.
 * Android Activity for Controlling Bluetooth LE Device Connectivity
 */

class DeviceControlActivity : Activity(), ActBle.ActBleListener {
    // Graphing Variables:
    private var mGraphInitializedBoolean = false
    private var mGraphAdapterCh1: GraphAdapter? = null
    private var mGraphAdapterCh2: GraphAdapter? = null
    private var mGraphAdapterCh3: GraphAdapter? = null
    private var mTimeDomainPlotAdapterCh1: XYPlotAdapter? = null
    private var mTimeDomainPlotAdapterCh2: XYPlotAdapter? = null
    private var mTimeDomainPlotAdapterCh3: XYPlotAdapter? = null
    private var mCh1: DataChannel? = null
    private var mCh2: DataChannel? = null
    private var mCh3: DataChannel? = null
    //Device Information
    private var mBleInitializedBoolean = false
    private lateinit var mBluetoothGattArray: Array<BluetoothGatt?>
    private var mLedWheelchairControlService: BluetoothGattService? = null
    //Classification
    private var mWheelchairControl = true //Default classifier.
    private var mActBle: ActBle? = null
    private var mDeviceName: String? = null
    private var mDeviceAddress: String? = null
    private var mConnected: Boolean = false
    private var mMSBFirst = false
    //Connecting to Multiple Devices
    private var deviceMacAddresses: Array<String>? = null
    private var mEEGConfigGattService: BluetoothGattService? = null
    private var mWheelchairGattIndex: Int = 0
    private var mEEGConfigGattIndex: Int = 0
    private var mEEGConnectedAllChannels = false
    // Classification
    private var mNumberPackets = -1
    private var mNumberOfClassifierCalls = 0
    private var mRunTrainingBool: Boolean = false
    //UI Elements - TextViews, Buttons, etc
    private var mBatteryLevel: TextView? = null
    private var mDataRate: TextView? = null
    private var mTrainingInstructions: TextView? = null
    private var mEMGClassText: TextView? = null
    private var mYfitTextView: TextView? = null
    private var mTogglePlots: ToggleButton? = null
    private var menu: Menu? = null
    //Data throughput counter
    private var mLastTime: Long = 0
    private var mLastTime2: Long = 0
    private var points = 0
    private var points2 = 0
    private val mTimerHandler = Handler()
    private var mTimerEnabled = false
    //Data Variables:
    private val batteryWarning = 20
    private var dataRate: Double = 0.toDouble()
    private var mScaleFactor = 20.0
    //Play Sound:
    private lateinit var mMediaBeep: MediaPlayer
    //Tensorflow:
    private var mTFRunModel = false
    private var mTFInferenceInterface: TensorFlowInferenceInterface? = null
    private var mTensorflowWindowSize = 128
    private val mClassificationBuffer = IntArray(5)

    //Drone Interface Stuff:
    private var mMiniDrone: MiniDrone? = null
    private var mDronePresent = false
    private var mDroneControl = false //Default classifier.
    private var mConnectionProgressDialog: ProgressDialog? = null
    private var mDownloadProgressDialog: ProgressDialog? = null
    private var mTakeOffLandBt: Button? = null
    private var mConnectDroneButton: Button? = null
    private var mDroneConnectionState = false
    private var mBatteryLevelDrone: TextView? = null
    private var mNbMaxDownload: Int = 0
    private var mCurrentDownloadIndex: Int = 0
    private var mARService: ARDiscoveryDeviceService? = null
    // Jumping Sumo Drone:
    private var mJSDrone: JSDrone? = null
    private val mJSDroneSpeedLR: Int = 20
    private val mJSDroneSpeedFWREV: Int = 20

    // EMG Training and Classification:
    private var mEMGClass = 0.0
    private val mEMGClassArray = DoubleArray(30000)
    private lateinit var mKNNParams: DoubleArray
    private var mClassifierReady = false
    private val yFitArray = DoubleArray(5)

    private val timeStamp: String
        get() = SimpleDateFormat("yyyy.MM.dd_HH.mm.ss", Locale.US).format(Date())

    private fun Double.format(digits: Int) = java.lang.String.format("%.${digits}f", this)!!

    private fun allValuesSame(y: IntArray): Boolean {
        return y.map { it == y[0] }.find { !it } == null
    }

    private fun lastThreeMatches(yfitarray: DoubleArray): Boolean {
        var b0 = false
        var b1 = false
        if (yfitarray[4] != 0.0) {
            b0 = yfitarray[4] == yfitarray[3]
            b1 = yfitarray[3] == yfitarray[2]
        }
        return b0 && b1
    }

    // Native Interface Function Handler:
    private val mNativeInterface = NativeInterfaceClass()

    //Bluetooth Classic - For Robotic Hand
    private var btAdapter: BluetoothAdapter? = null
    private var btSocket: BluetoothSocket? = null

    private lateinit var mConnectedThread: ConnectedThread
    private val sb = StringBuilder()
    private var flag = 0
    //Bluetooth Classic - For Robotic Hand
//    var mHandler = Handler()

    public override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activit_dev_ctrl_alt)
        //Set orientation of device based on screen type/size:
        requestedOrientation = ActivityInfo.SCREEN_ORIENTATION_LANDSCAPE
        //Receive Intents:
        //Receive Intents:
        val intent = intent
        deviceMacAddresses = intent.getStringArrayExtra(MainActivity.INTENT_DEVICES_KEY)
        val deviceDisplayNames = intent.getStringArrayExtra(MainActivity.INTENT_DEVICES_NAMES)
        val stimulusSeconds = intent.getStringArrayExtra(MainActivity.INTENT_DELAY_VALUE_SECONDS)
        if (intent.extras != null) {
            mRunTrainingBool = intent.extras!!.getBoolean(MainActivity.INTENT_TRAIN_BOOLEAN)
            Log.e(TAG, "Train Activity (true/false): $mRunTrainingBool")
        } else {
            Log.e(TAG, "ERROR: intent.getExtras = null")
        }
        mScaleFactor = Integer.valueOf(stimulusSeconds[0])!!.toDouble()
        mDeviceName = deviceDisplayNames[0]
        mDeviceAddress = deviceMacAddresses!![0]
        Log.d(TAG, "Device Names: " + Arrays.toString(deviceDisplayNames))
        Log.d(TAG, "Device MAC Addresses: " + Arrays.toString(deviceMacAddresses))
        Log.d(TAG, Arrays.toString(deviceMacAddresses))
        //Set up action bar:
        if (actionBar != null) {
            actionBar!!.setDisplayHomeAsUpEnabled(true)
        }
        val actionBar = actionBar
        actionBar!!.setBackgroundDrawable(ColorDrawable(Color.parseColor("#6078ef")))
        //Flag to keep screen on (stay-awake):
        window.addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON)
        //Set up TextViews
        val mExportButton = findViewById<Button>(R.id.button_export)
        mConnectDroneButton = findViewById(R.id.droneConnectButton)
        mBatteryLevel = findViewById(R.id.batteryText)
        mDataRate = findViewById(R.id.dataRate)
        mTrainingInstructions = findViewById(R.id.trainingInstructions)
        mEMGClassText = findViewById(R.id.emgClassText)
        mYfitTextView = findViewById(R.id.textViewYfit)
        mTakeOffLandBt = findViewById(R.id.buttonS)
        mDataRate!!.text = "..."
        val ab = getActionBar()
        ab!!.title = mDeviceName
        ab.subtitle = mDeviceAddress
        //Initialize Bluetooth
        if (!mBleInitializedBoolean) initializeBluetoothArray()
        mMediaBeep = MediaPlayer.create(this, R.raw.beep_01a)
        mLastTime = System.currentTimeMillis()
        mLastTime2 = System.currentTimeMillis()
        updateTrainingView(mRunTrainingBool)
        mTogglePlots = findViewById(R.id.toggleButtonCh1)
        mTogglePlots!!.setOnCheckedChangeListener { _, b ->
            if (!b) {
                mGraphAdapterCh3?.clearPlot()
                mGraphAdapterCh2?.clearPlot()
                mGraphAdapterCh1?.clearPlot()
            }
            mGraphAdapterCh1!!.plotData = b
            mGraphAdapterCh2!!.plotData = b
            mGraphAdapterCh3!!.plotData = b
        }
        resetScaleButton.setOnClickListener {
            if (mZAccValue < 0.70) {
                mZAccMaxThreshold = mZAccValue /*current value*/ + 0.20
                mZAccMinThreshold = mZAccValue - 0.20
                val s = "Max = ${mZAccMaxThreshold.format(2)} \n" +
                        " Min = ${mZAccMinThreshold.format(2)}"
                Toast.makeText(applicationContext, s, Toast.LENGTH_SHORT).show()
            } else {
                Toast.makeText(applicationContext, "No changes made! \n Check Your position", Toast.LENGTH_LONG).show()
            }
        }
        val toggleButton1 = findViewById<ToggleButton>(R.id.toggleButtonDroneControl)
        toggleButton1.setOnCheckedChangeListener { _, b ->
            mDroneControl = b
        }
        mExportButton.setOnClickListener { exportData() }
        findViewById<View>(R.id.buttonFwd).setOnTouchListener { v, event ->
            when (event.action) {
                MotionEvent.ACTION_DOWN -> {
                    v.isPressed = true
                    sendDroneCommand(1)
                    executeWheelchairCommand(1)
                }

                MotionEvent.ACTION_UP -> {
                    mMiniDrone?.setPitch(0.toByte())
                    mMiniDrone?.setFlag(0.toByte())
                    sendDroneCommand(0)
                    executeWheelchairCommand(0)
                }
                else -> {
                }
            }
            true
        }
        findViewById<View>(R.id.buttonR).setOnTouchListener { v, event ->
            when (event.action) {
                MotionEvent.ACTION_DOWN -> {
                    v.isPressed = true
//                    mMiniDrone?.setYaw(50.toByte())
                    sendDroneCommand(2)
                    executeWheelchairCommand(3)
                    processClassifiedData(0.0)
                }

                MotionEvent.ACTION_UP -> {
                    v.isPressed = false
                    mMiniDrone?.setYaw(0.toByte())
                    sendDroneCommand(0)
                    executeWheelchairCommand(0)
                }
                else -> {
                }
            }
            true
        }
        buttonL.setOnTouchListener { v, event ->
            when (event.action) {
                MotionEvent.ACTION_DOWN -> {
                    v.isPressed = true
                    sendDroneCommand(3)
                    processClassifiedData(1.0)
                }

                MotionEvent.ACTION_UP -> {
                    v.isPressed = false
                    sendDroneCommand(0)
                }
                else -> {
                }
            }
            true
        }
        mTakeOffLandBt?.setOnClickListener {
            if (mMiniDrone != null) {
                when (mMiniDrone?.flyingState) {
                    ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_LANDED -> mMiniDrone?.takeOff()
                    ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_FLYING, ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_HOVERING -> mMiniDrone?.land()
                    else -> {/*Do Nothing*/
                    }
                }
            }
        }
        mConnectDroneButton?.setOnClickListener {
            if (!mDroneConnectionState) {
                if (connectDrone()) {
                    mDroneConnectionState = true
                    val s = "Disconnect"
                    mConnectDroneButton?.text = s
                }
            } else {
                if (disconnectDrone()) {
                    mDroneConnectionState = false
                    val s = "Connect"
                    mConnectDroneButton?.text = s
                }
            }
            Log.d(TAG, "onClick: mDroneConnectionState: $mDroneConnectionState")
        }
        mARService = intent.getParcelableExtra(MainActivity.EXTRA_DRONE_SERVICE)
        //Attempt to connect to BT Hand Automatically
//        btAdapter = BluetoothAdapter.getDefaultAdapter()       // get Bluetooth adapter
//        checkBTState()
//        mHandler = object : Handler() {
//            override fun handleMessage(msg: android.os.Message) {
//                when (msg.what) {
//                    RECIEVE_MESSAGE -> {
//                        val readBuf = msg.obj as ByteArray
//                        val strIncom = String(readBuf, 0, msg.arg1)
//                        sb.append(strIncom)
//                        val endOfLineIndex = sb.indexOf("\r\n")
//                        if (endOfLineIndex > 0) {
//                            val sbprint = sb.substring(0, endOfLineIndex)
//                            sb.delete(0, sb.length)
//                            flag++
//                            Log.i(TAG, "flag: " + flag.toString() + "Sbprint: " + sbprint)
//                        }
//                    }
//                }
//            }
//        }
//        connectToClassicBTHand()
    }

    private fun checkBTState() {
        // Check for Bluetooth support and then check to make sure it is turned on
        // Emulator doesn't support Bluetooth and will return null
        if (btAdapter == null) {
            Log.e(TAG, "Fatal Error: " + "Bluetooth not support")
        } else {
            if (btAdapter!!.isEnabled) {
                Log.d(TAG, "...Bluetooth ON...")
            } else {
                //Prompt user to turn on Bluetooth
                val enableBtIntent = Intent(BluetoothAdapter.ACTION_REQUEST_ENABLE)
                startActivityForResult(enableBtIntent, 1)
            }
        }
    }

    @Throws(IOException::class)
    private fun createBluetoothSocket(device: BluetoothDevice): BluetoothSocket {
        try {
            val m = device.javaClass.getMethod("createInsecureRfcommSocketToServiceRecord", UUID::class.java)
            return m.invoke(device, MY_UUID) as BluetoothSocket
        } catch (e: Exception) {
            Log.e(TAG, "Could not create Insecure RFComm Connection", e)
        }

        return device.createRfcommSocketToServiceRecord(MY_UUID)
    }

    private fun connectToClassicBTHand() {
        val device = btAdapter!!.getRemoteDevice(address)
        try {
            btSocket = createBluetoothSocket(device)
        } catch (e: IOException) {
            Log.e(TAG, "socketCreate fail" + e.message)
        }

        btAdapter!!.cancelDiscovery()
        Log.d(TAG, "...Connecting...")
        try {
            btSocket!!.connect()
            Log.d(TAG, "....Connection ok...")
        } catch (e: IOException) {
            try {
                btSocket!!.close()
            } catch (e2: IOException) {
                Log.e(TAG, "In onResume() and unable to close socket during connection failure" + e2.message + ".")
            }

        }

        // Create a data stream so we can talk to server.
        Log.d(TAG, "...Create Socket...")

        mConnectedThread = ConnectedThread(btSocket)
        mConnectedThread.start()
    }

    private fun disconnectDrone(): Boolean {
        if (mMiniDrone != null) {
            mConnectionProgressDialog = ProgressDialog(this, R.style.AppCompatAlertDialogStyle)
            mConnectionProgressDialog?.isIndeterminate = true
            mConnectionProgressDialog?.setMessage("Disconnecting ...")
            mConnectionProgressDialog?.setCancelable(false)
            mConnectionProgressDialog?.show()
            if (!mMiniDrone!!.disconnect()) {
                finish()
                return false
            }
        }
        return true
    }

    private fun connectDrone(): Boolean {
        if (mARService != null) {
            mBatteryLevelDrone?.visibility = View.VISIBLE
            mDronePresent = true
            // TODO: How to check drone type:
            val droneName = mARService?.name?.toLowerCase()
            Log.e(TAG, "mARService->name: $droneName")
            if (droneName!!.contains("mambo")) {
                mMiniDrone = MiniDrone(this, mARService!!)
                mMiniDrone?.addListener(mMiniDroneListener)
                if (mMiniDrone != null && ARCONTROLLER_DEVICE_STATE_ENUM.ARCONTROLLER_DEVICE_STATE_RUNNING != mMiniDrone?.getConnectionState()) {
                    mConnectionProgressDialog = ProgressDialog(this, R.style.AppCompatAlertDialogStyle)
                    mConnectionProgressDialog?.isIndeterminate = true
                    mConnectionProgressDialog?.setMessage("Connecting ...")
                    mConnectionProgressDialog?.setCancelable(false)
                    mConnectionProgressDialog?.show()
                    // if the connection to the MiniDrone fails, finish the activity
                    if (!mMiniDrone!!.connect()) {
                        Toast.makeText(applicationContext, "Failed to Connect Drone", Toast.LENGTH_LONG).show()
                        finish()
                        return false
                    }
                    return true
                } else {
                    return false
                }
            } else if (droneName.contains("diesel") || droneName.contains("sumo") || droneName.contains("black")) {
                buttonS.visibility = View.GONE
                mJSDrone = JSDrone(this, mARService!!)
                mJSDrone?.addListener(mJSListener)
                if (ARCONTROLLER_DEVICE_STATE_ENUM.ARCONTROLLER_DEVICE_STATE_RUNNING != mJSDrone!!.connectionState) {
                    mConnectionProgressDialog = ProgressDialog(this, R.style.AppCompatAlertDialogStyle)
                    mConnectionProgressDialog!!.isIndeterminate = true
                    mConnectionProgressDialog!!.setMessage("Connecting ...")
                    mConnectionProgressDialog!!.show()

                    // if the connection to the Jumping fails, finish the activity
                    if (!mJSDrone!!.connect()) {
                        Toast.makeText(applicationContext, "Failed to Connect Drone", Toast.LENGTH_LONG).show()
                        finish()
                        return false
                    }
                    return true
                } else {
                    return false
                }
            }
        } else {
            Log.e(TAG, "Drone Service is Null")
            return false
        }
        return false
    }

    private fun executeWheelchairCommand(command: Int) {
        val bytes = ByteArray(1)
        when (command) {
            0 -> bytes[0] = 0x00.toByte()
            1 -> bytes[0] = 0x01.toByte() // Stop
            2 -> bytes[0] = 0xF0.toByte() // Rotate Left
            3 -> bytes[0] = 0x0F.toByte() // Rotate Right ??
            4 -> bytes[0] = 0xFF.toByte() // TODO: 6/27/2017 Disconnect instead of reverse?
            else -> {
            }
        }
        if (mLedWheelchairControlService != null && mWheelchairControl) {
            Log.e(TAG, "SendingCommand: $command")
            Log.e(TAG, "SendingCommand (byte): " + DataChannel.byteArrayToHexString(bytes))
            mActBle!!.writeCharacteristic(mBluetoothGattArray[mWheelchairGattIndex]!!, mLedWheelchairControlService!!.getCharacteristic(AppConstant.CHAR_WHEELCHAIR_CONTROL), bytes)
        }
    }

    private fun sendDroneCommand(command: Int) {
        if (mDronePresent && mDroneControl && mMiniDrone != null) {
            Log.e(TAG, "SendingCommand to MiniDrone: $command")
            when (command) {
                0 -> {
                    // Do nothing/Hover //Reset conditions:
                    mMiniDrone?.setYaw(0.toByte())
                    mMiniDrone?.setPitch(0.toByte())
                    mMiniDrone?.setFlag(0.toByte())
                }
                1 -> {
                    // Move Fwd:
                    mMiniDrone?.setPitch(50.toByte())
                    mMiniDrone?.setFlag(1.toByte())
                }
                2 -> //RotR
                    mMiniDrone?.setYaw(50.toByte())
                3 -> // Take off or Land
                    when (mMiniDrone?.flyingState) {
                        ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_LANDED -> mMiniDrone?.takeOff()
                        ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_FLYING,
                        ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_HOVERING -> mMiniDrone?.land()
                        else -> {
                        } //Do nothing.
                    }
                4 -> { // Roll Left
                    mMiniDrone?.setRoll((-50).toByte())
                    mMiniDrone?.setFlag(1.toByte())
                }
                5 -> { // Roll Right
                    mMiniDrone?.setRoll(50.toByte())
                    mMiniDrone?.setFlag(1.toByte())
                }
                else -> {
                }
            }
        }
        if (mDronePresent && mDroneControl && mJSDrone != null) {
            Log.e(TAG, "SendingCommand to JSDrone: $command")
            when (command) {
                0 -> { //Do Nothing
                    mJSDrone?.setTurn(0.toByte())
                    mJSDrone?.setSpeed(0.toByte())
                    mJSDrone?.setFlag(0.toByte())
                }
                1 -> {  //FWD:
                    mJSDrone?.setSpeed(mJSDroneSpeedFWREV.toByte())
                    mJSDrone?.setTurn(0.toByte())
                    mJSDrone?.setFlag(1.toByte())
                }
                2 -> { // RIGHT
                    mJSDrone?.setSpeed(0.toByte())
                    mJSDrone?.setTurn((mJSDroneSpeedLR).toByte())
                    mJSDrone?.setFlag(1.toByte())
                }
                3 -> { // LEFT
                    mJSDrone?.setSpeed(0.toByte())
                    mJSDrone?.setTurn((-mJSDroneSpeedLR).toByte())
                    mJSDrone?.setFlag(1.toByte())
                }
                4, 5 -> {
                    // No Change in State (ignore command).
                }
                else -> {
                    mJSDrone?.setTurn(0.toByte())
                    mJSDrone?.setSpeed(0.toByte())
                    mJSDrone?.setFlag(0.toByte())
                }
            }
        }
    }

    private fun exportData() {
        try {
            terminateDataFileWriter()
        } catch (e: IOException) {
            Log.e(TAG, "IOException in saveDataFile")
            e.printStackTrace()
        }
        val files = ArrayList<Uri>()
        val context = applicationContext
        val uii = FileProvider.getUriForFile(context, context.packageName + ".provider", mPrimarySaveDataFile!!.file)
        files.add(uii)
        if (mSaveFileMPU != null) {
            val uii2 = FileProvider.getUriForFile(context, context.packageName + ".provider", mSaveFileMPU!!.file)
            files.add(uii2)
        }
        val exportData = Intent(Intent.ACTION_SEND_MULTIPLE)
        exportData.putExtra(Intent.EXTRA_SUBJECT, "EMG/MPU Sensor Data Export Details")
        exportData.putParcelableArrayListExtra(Intent.EXTRA_STREAM, files)
        exportData.type = "text/html"
        startActivity(exportData)
    }


    @Throws(IOException::class)
    private fun terminateDataFileWriter() {
        mPrimarySaveDataFile?.terminateDataFileWriter()
    }

    public override fun onResume() {
        mNativeInterface.jmainInitialization(false)
        if (mRedrawer != null) {
            mRedrawer!!.start()
        }
        super.onResume()
    }

    override fun onPause() {
        if (mRedrawer != null) mRedrawer!!.pause()
        super.onPause()
    }

    private fun initializeBluetoothArray() {
        val mBluetoothManager = getSystemService(Context.BLUETOOTH_SERVICE) as BluetoothManager
        val mBluetoothDeviceArray = arrayOfNulls<BluetoothDevice>(deviceMacAddresses!!.size)
        Log.d(TAG, "Device Addresses: " + Arrays.toString(deviceMacAddresses))
        if (deviceMacAddresses != null) {
            for (i in deviceMacAddresses!!.indices) {
                mBluetoothDeviceArray[i] = mBluetoothManager.adapter.getRemoteDevice(deviceMacAddresses!![i])
            }
        } else {
            Log.e(TAG, "No Devices Queued, Restart!")
            Toast.makeText(this, "No Devices Queued, Restart!", Toast.LENGTH_SHORT).show()
        }
        mActBle = ActBle(this, mBluetoothManager, this)
        mBluetoothGattArray = Array(deviceMacAddresses!!.size) { i -> mActBle!!.connect(mBluetoothDeviceArray[i]) }
        for (i in mBluetoothDeviceArray.indices) {
            Log.e(TAG, "Connecting to Device: Name: " + (mBluetoothDeviceArray[i]!!.name + " \nMAC:" + mBluetoothDeviceArray[i]!!.address))
            if ("WheelchairControl" == mBluetoothDeviceArray[i]!!.name) {
                mWheelchairGattIndex = i
                Log.e(TAG, "mWheelchairGattIndex: $mWheelchairGattIndex")
                continue //we are done initializing
            } else {
                mEEGConfigGattIndex = i
            }

            val btDeviceName = mBluetoothDeviceArray[i]?.name?.toLowerCase()
            mMSBFirst = when {
                btDeviceName == null -> false
                btDeviceName.contains("EMG 250Hz") -> false
                btDeviceName.contains("nrf52") -> true
                btDeviceName.contains("emg_") -> true
                else -> false
            }
            mSampleRate = when {
                btDeviceName == null -> 250
                btDeviceName.contains("8k") -> 8000
                btDeviceName.contains("4k") -> 4000
                btDeviceName.contains("2k") -> 2000
                btDeviceName.contains("1k") -> 1000
                btDeviceName.contains("500") -> 500
                else -> 250
            }
//            mPacketBuffer = mSampleRate / 250
            Log.e(TAG, "mSampleRate: " + mSampleRate + "Hz")
            if (!mGraphInitializedBoolean) setupGraph()

            mGraphAdapterCh1!!.setxAxisIncrementFromSampleRate(mSampleRate)
            mGraphAdapterCh2!!.setxAxisIncrementFromSampleRate(mSampleRate)
            mGraphAdapterCh3!!.setxAxisIncrementFromSampleRate(mSampleRate)

            mGraphAdapterCh1!!.setSeriesHistoryDataPoints(250 * 5)
            mGraphAdapterCh2!!.setSeriesHistoryDataPoints(250 * 5)
            mGraphAdapterCh3!!.setSeriesHistoryDataPoints(250 * 5)
            val fileNameTimeStamped = "EMG_3ChData_" + timeStamp + "_" + mSampleRate.toString() + "Hz"
            Log.e(TAG, "fileTimeStamp: $fileNameTimeStamped")
            try {
                mPrimarySaveDataFile = SaveDataFile("/EMGData", fileNameTimeStamped,
                        24, 1.toDouble() / mSampleRate)
            } catch (e: IOException) {
                Log.e(TAG, "initializeBluetoothArray: IOException", e)
            }

            mPrimarySaveDataFile!!.setSaveTimestamps(false)
            mPrimarySaveDataFile!!.setFpPrecision(64)
            mPrimarySaveDataFile!!.setIncludeClass(true)
        }
        mBleInitializedBoolean = true
    }

    private fun createNewFileMPU() {
        val directory = "/MPUData"
        val fileNameTimeStamped = "MPUData_$timeStamp"
        if (mSaveFileMPU == null) {
            Log.e(TAG, "fileTimeStamp: $fileNameTimeStamped")
            mSaveFileMPU = SaveDataFile(directory, fileNameTimeStamped,
                    16, 0.020)
        } else if (!mSaveFileMPU!!.initialized) {
//            Log.e(TAG, "New Filename: " + fileNameTimeStamped)
//            mSaveFileMPU?.createNewFile(directory, fileNameTimeStamped)
        }
    }


    private fun setupGraph() {
        // Initialize our XYPlot reference:
        mGraphAdapterCh1 = GraphAdapter(mSampleRate * 4, "EMG Data Ch 1", false, Color.BLUE)
        mGraphAdapterCh2 = GraphAdapter(mSampleRate * 4, "EMG Data Ch 2", false, Color.RED)
        mGraphAdapterCh3 = GraphAdapter(mSampleRate * 4, "EMG Data Ch 3", false, Color.GRAY)
        //PLOT By default
        mGraphAdapterCh1!!.plotData = true
        mGraphAdapterCh2!!.plotData = true
        mGraphAdapterCh3!!.plotData = true
        mGraphAdapterCh1!!.setPointWidth(2.toFloat())
        mGraphAdapterCh2!!.setPointWidth(2.toFloat())
        mGraphAdapterCh3!!.setPointWidth(2.toFloat())

        mTimeDomainPlotAdapterCh1 = XYPlotAdapter(findViewById(R.id.emgTimeDomainXYPlot), false, 1000)
        if (mTimeDomainPlotAdapterCh1!!.xyPlot != null) {
            mTimeDomainPlotAdapterCh1!!.xyPlot!!.addSeries(mGraphAdapterCh1!!.series, mGraphAdapterCh1!!.lineAndPointFormatter)
        }
        mTimeDomainPlotAdapterCh2 = XYPlotAdapter(findViewById(R.id.emgTimeDomainXYPlot2), false, 1000)
        if (mTimeDomainPlotAdapterCh2!!.xyPlot != null) {
            mTimeDomainPlotAdapterCh2!!.xyPlot!!.addSeries(mGraphAdapterCh2!!.series, mGraphAdapterCh2!!.lineAndPointFormatter)
        }
        mTimeDomainPlotAdapterCh3 = XYPlotAdapter(findViewById(R.id.emgTimeDomainXYPlot3), false, 1000)
        if (mTimeDomainPlotAdapterCh3!!.xyPlot != null) {
            mTimeDomainPlotAdapterCh3?.xyPlot?.addSeries(mGraphAdapterCh3?.series, mGraphAdapterCh3?.lineAndPointFormatter)
        }
        val xyPlotList = listOf(mTimeDomainPlotAdapterCh1!!.xyPlot, mTimeDomainPlotAdapterCh2!!.xyPlot, mTimeDomainPlotAdapterCh3?.xyPlot)
        mRedrawer = Redrawer(xyPlotList, 30f, false)
        mRedrawer!!.start()
        mGraphInitializedBoolean = true
    }

    private fun setNameAddress(name_action: String?, address_action: String?) {
        val name = menu!!.findItem(R.id.action_title)
        val address = menu!!.findItem(R.id.action_address)
        name.title = name_action
        address.title = address_action
        invalidateOptionsMenu()
    }

    override fun onDestroy() {
        mRedrawer?.finish()
        disconnectAllBLE()
        try {
            terminateDataFileWriter()
        } catch (e: IOException) {
            Log.e(TAG, "IOException in saveDataFile")
            e.printStackTrace()
        }
        if (mJSDrone != null) {
            mJSDrone?.dispose()
        }
        if (mMiniDrone != null) {
            mMiniDrone?.dispose()
        }
        stopMonitoringRssiValue()
        mNativeInterface.jmainInitialization(true) //Just a technicality, doesn't actually do anything
        super.onDestroy()
        if (btSocket != null) {
            try {
                btSocket!!.close()
            } catch (e2: IOException) {
                Log.e(TAG, "In onPause() and failed to close socket." + e2.message + ".")
            }
        }
    }

    private fun disconnectAllBLE() {
        if (mActBle != null) {
            for (bluetoothGatt in mBluetoothGattArray) {
                mActBle!!.disconnect(bluetoothGatt!!)
                mConnected = false
                resetMenuBar()
            }
        }
    }

    private fun resetMenuBar() {
        runOnUiThread {
            if (menu != null) {
                menu!!.findItem(R.id.menu_connect).isVisible = true
                menu!!.findItem(R.id.menu_disconnect).isVisible = false
            }
        }
    }

    override fun onCreateOptionsMenu(menu: Menu): Boolean {
        menuInflater.inflate(R.menu.menu_device_control, menu)
        menuInflater.inflate(R.menu.actionbar_item, menu)
        if (mConnected) {
            menu.findItem(R.id.menu_connect).isVisible = false
            menu.findItem(R.id.menu_disconnect).isVisible = true
        } else {
            menu.findItem(R.id.menu_connect).isVisible = true
            menu.findItem(R.id.menu_disconnect).isVisible = false
        }
        this.menu = menu
        setNameAddress(mDeviceName, mDeviceAddress)
        return true
    }

    override fun onOptionsItemSelected(item: MenuItem): Boolean {
        when (item.itemId) {
            R.id.menu_connect -> {
                if (mActBle != null) {
                    initializeBluetoothArray()
                }
                connect()
                return true
            }
            R.id.menu_disconnect -> {
                if (mActBle != null) {
                    disconnectAllBLE()
                }
                return true
            }
            android.R.id.home -> {
                if (mActBle != null) {
                    disconnectAllBLE()
                }
                NavUtils.navigateUpFromSameTask(this)
                onBackPressed()
                return true
            }
            R.id.action_settings -> {
                launchSettingsMenu()
                return true
            }
            R.id.action_export -> {
                exportData()
                return true
            }
        }
        return super.onOptionsItemSelected(item)
    }

    override fun onActivityResult(requestCode: Int, resultCode: Int, data: Intent?) {
        if (requestCode == 1) {
            val context = applicationContext
            //File Save Stuff
            val saveTimestamps = PreferencesFragment.saveTimestamps(context)
            val precision = (if (PreferencesFragment.setBitPrecision(context)) 64 else 32).toShort()
            val saveClass = PreferencesFragment.saveClass(context)
            mPrimarySaveDataFile!!.setSaveTimestamps(saveTimestamps)
            mPrimarySaveDataFile!!.setFpPrecision(precision)
            mPrimarySaveDataFile!!.setIncludeClass(saveClass)
            val filterData = PreferencesFragment.setFilterData(context)
            //TODO: for now just ch1:
            if (mGraphAdapterCh1 != null) {
                mFilterData = filterData
            }
            /**
             * Settings for ADS1299:
             */
            val registerConfigBytes = Arrays.copyOf(ADS1299_DEFAULT_BYTE_CONFIG, ADS1299_DEFAULT_BYTE_CONFIG.size)
            when (PreferencesFragment.setSampleRate(context)) {
                0 -> {
                    registerConfigBytes[0] = 0x96.toByte()
                }
                1 -> {
                    registerConfigBytes[0] = 0x95.toByte()
                }
                2 -> {
                    registerConfigBytes[0] = 0x94.toByte()
                }
                3 -> {
                    registerConfigBytes[0] = 0x93.toByte()
                }
                4 -> {
                    registerConfigBytes[0] = 0x92.toByte()
                }
            }
            val numChEnabled = PreferencesFragment.setNumberChannelsEnabled(context)
            Log.e(TAG, "numChEnabled: " + numChEnabled.toString())
            // Set
            when (numChEnabled) {
                1 -> {
                    registerConfigBytes[12] = 0b0000_0001
                    registerConfigBytes[13] = 0b0000_0001
                }
                2 -> {
                    registerConfigBytes[12] = 0b0000_0011
                    registerConfigBytes[13] = 0b0000_0011
                }
                3 -> {
                    registerConfigBytes[12] = 0b0000_0111
                    registerConfigBytes[13] = 0b0000_0111
                }
                4 -> {
                    registerConfigBytes[12] = 0b0000_1111
                    registerConfigBytes[13] = 0b0000_1111
                }
            }

            //Set all to disable.
            for (i in 4..7) registerConfigBytes[i] = 0xE1.toByte()
            Log.e(TAG, "SettingsNew0: " + DataChannel.byteArrayToHexString(registerConfigBytes))
            for (i in 4..(3 + numChEnabled)) {
                registerConfigBytes[i] = 0x00.toByte()
            }
            Log.e(TAG, "SettingsNew1: " + DataChannel.byteArrayToHexString(registerConfigBytes))
            val gain12 = PreferencesFragment.setGainCh12(context) //Check if ch enabled first
            for (i in 4..5) { //Checks first bit enabled on chs 1 & 2.
                registerConfigBytes[i] = when (registerConfigBytes[i] and 0x80.toByte()) {
                    0x80.toByte() -> registerConfigBytes[i] // do nothing if ch disabled
                    else -> registerConfigBytes[i] or (gain12 shl 4).toByte() // Shift gain into 0xxx_0000 position
                }
            }
            if (PreferencesFragment.setSRB1(context)) {
                registerConfigBytes[20] = 0x20.toByte()
                registerConfigBytes[13] = 0b0000_0000 // Turn off BIAS_SENSN with SRB 1
            } else {
                registerConfigBytes[20] = 0x00.toByte()
            }
            Log.e(TAG, "SettingsNew: " + DataChannel.byteArrayToHexString(registerConfigBytes))
            writeNewADS1299Settings(registerConfigBytes)
        }
        super.onActivityResult(requestCode, resultCode, data)
    }

    private fun writeNewADS1299Settings(bytes: ByteArray) {
        Log.e(TAG, "bytesOriginal: " + DataChannel.byteArrayToHexString(ADS1299_DEFAULT_BYTE_CONFIG))
        if (mEEGConfigGattService != null) {
            Log.e(TAG, "SendingCommand (byte): " + DataChannel.byteArrayToHexString(bytes))
            mActBle!!.writeCharacteristic(mBluetoothGattArray[mEEGConfigGattIndex]!!, mEEGConfigGattService!!.getCharacteristic(AppConstant.CHAR_EEG_CONFIG), bytes)
            //Should notify/update after writing
        }
    }

    private fun launchSettingsMenu() {
        val intent = Intent(applicationContext, SettingsActivity::class.java)
        startActivityForResult(intent, 1)
    }

    private fun connect() {
        runOnUiThread {
            val menuItem = menu!!.findItem(R.id.action_status)
            menuItem.title = "Connecting..."
        }
    }

    override fun onServicesDiscovered(gatt: BluetoothGatt, status: Int) {
        Log.i(TAG, "onServicesDiscovered")
        if (status == BluetoothGatt.GATT_SUCCESS) {
            for (service in gatt.services) {
                if (service == null || service.uuid == null) {
                    continue
                }

                if (AppConstant.SERVICE_WHEELCHAIR_CONTROL == service.uuid) {
                    mLedWheelchairControlService = service
                    Log.i(TAG, "BLE Wheelchair Control Service found")
                }

                if (AppConstant.SERVICE_DEVICE_INFO == service.uuid) {
                    //Read the device serial number (if available)
                    if (service.getCharacteristic(AppConstant.CHAR_SERIAL_NUMBER) != null) {
                        mActBle!!.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_SERIAL_NUMBER))
                    }
                    //Read the device software version (if available)
                    if (service.getCharacteristic(AppConstant.CHAR_SOFTWARE_REV) != null) {
                        mActBle!!.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_SOFTWARE_REV))
                    }
                }

                if (AppConstant.SERVICE_EEG_SIGNAL == service.uuid) {
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CONFIG) != null) {
                        mEEGConfigGattService = service
                        mActBle!!.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CONFIG))
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CONFIG), true)
                    }

                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH1_SIGNAL) != null) {
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH1_SIGNAL), true)
                        if (mCh1 == null) mCh1 = DataChannel(false, mMSBFirst, 120 * mSampleRate)
                    }
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH2_SIGNAL) != null) {
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH2_SIGNAL), true)
                        if (mCh2 == null) mCh2 = DataChannel(false, mMSBFirst, 120 * mSampleRate)
                    }
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH3_SIGNAL) != null) {
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH3_SIGNAL), true)
                        if (mCh3 == null) mCh3 = DataChannel(false, mMSBFirst, 120 * mSampleRate)
                    }
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH4_SIGNAL) != null)
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH4_SIGNAL), true)
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH5_SIGNAL) != null)
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH5_SIGNAL), true)
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH6_SIGNAL) != null)
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH6_SIGNAL), true)
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH7_SIGNAL) != null)
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH7_SIGNAL), true)
                    if (service.getCharacteristic(AppConstant.CHAR_EEG_CH8_SIGNAL) != null)
                        mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH8_SIGNAL), true)
                }

                if (AppConstant.SERVICE_BATTERY_LEVEL == service.uuid) { //Read the device battery percentage
                    mActBle!!.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_BATTERY_LEVEL))
                    mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_BATTERY_LEVEL), true)
                }

                if (AppConstant.SERVICE_MPU == service.uuid) {
                    mActBle!!.setCharacteristicNotifications(gatt, service.getCharacteristic(AppConstant.CHAR_MPU_COMBINED), true)
                    mMPU = DataChannel(false, true, 0)
                    mAccX = DataChannel(false, true, 50)
                    mAccY = DataChannel(false, true, 50)
                    mAccZ = DataChannel(false, true, 50)
                    createNewFileMPU()
                }
            }
            //Run process only once:
            mActBle?.runProcess()
        }
    }

    override fun onCharacteristicRead(gatt: BluetoothGatt, characteristic: BluetoothGattCharacteristic, status: Int) {
        Log.i(TAG, "onCharacteristicRead")
        if (status == BluetoothGatt.GATT_SUCCESS) {
            if (AppConstant.CHAR_BATTERY_LEVEL == characteristic.uuid) {
                if (characteristic.value != null) {
                    val batteryLevel = characteristic.getIntValue(BluetoothGattCharacteristic.FORMAT_UINT16, 0)
                    updateBatteryStatus(batteryLevel)
                    Log.i(TAG, "Battery Level :: $batteryLevel")
                }
            }
            //TODO: NEED TO CHANGE mSampleRate, DataChannel[], and GraphAdapter[] here.
            if (AppConstant.CHAR_EEG_CONFIG == characteristic.uuid) {
                if (characteristic.value != null) {
                    val readValue = characteristic.value
                    Log.e(TAG, "onCharacteriticRead: \n" +
                            "CHAR_EEG_CONFIG: " + DataChannel.byteArrayToHexString(readValue))
                }
            }
        } else {
            Log.e(TAG, "onCharacteristic Read Error$status")
        }
    }

    override fun onCharacteristicChanged(gatt: BluetoothGatt, characteristic: BluetoothGattCharacteristic) {
        if (AppConstant.CHAR_EEG_CONFIG == characteristic.uuid) {
            if (characteristic.value != null) {
                val readValue = characteristic.value
                Log.e(TAG, "onCharacteristicChanged: \n" +
                        "CHAR_EEG_CONFIG: " + DataChannel.byteArrayToHexString(readValue))
                when (readValue[0] and 0x0F.toByte()) {
                    0x06.toByte() -> mSampleRate = 250
                    0x05.toByte() -> mSampleRate = 500
                    0x04.toByte() -> mSampleRate = 1000
                    0x03.toByte() -> mSampleRate = 2000
                    0x02.toByte() -> mSampleRate = 4000
                }
                //RESET mCH1 & mCH2:
                mCh1?.classificationBufferSize = 4 * mSampleRate
                Log.e(TAG, "Updated Sample Rate: $mSampleRate")
            }
        }

        if (AppConstant.CHAR_BATTERY_LEVEL == characteristic.uuid) {
            val batteryLevel = characteristic.getIntValue(BluetoothGattCharacteristic.FORMAT_UINT16, 0)!!
            updateBatteryStatus(batteryLevel)
        }

        /**
         * Either Channel 1 or Channel 2 Enabled (depending on circuit config).
         */
        if (AppConstant.CHAR_EEG_CH1_SIGNAL == characteristic.uuid) {
            val mNewEEGdataBytes = characteristic.value
            if (!mCh1!!.chEnabled) {
                mCh1!!.chEnabled = true
            }
            getDataRateBytes(mNewEEGdataBytes.size)
            if (mEEGConnectedAllChannels) {
                mCh1!!.handleNewData(mNewEEGdataBytes, false)
                if (mCh1!!.packetCounter.toInt() == mPacketBuffer) {
                    addToGraphBuffer(mCh1!!, mGraphAdapterCh1, true)
                    //TODO: Update Training Routine
//                    if (mNumberPackets % 10 == 0 && mClassifierReady) {
//                        val mClassifyThread = Thread(mClassifyEMGThread)
//                        mClassifyThread.start()
//                    }
                }
            }
        }

        if (AppConstant.CHAR_EEG_CH2_SIGNAL == characteristic.uuid) {
            if (!mCh2!!.chEnabled) {
                mCh2!!.chEnabled = true
            }
            val mNewEEGdataBytes = characteristic.value
            val byteLength = mNewEEGdataBytes.size
            getDataRateBytes(byteLength)
            if (mEEGConnectedAllChannels) {
                mCh2!!.handleNewData(mNewEEGdataBytes, false)
                if (mCh2!!.packetCounter.toInt() == mPacketBuffer) {
                    addToGraphBuffer(mCh2!!, mGraphAdapterCh2, false)
                }
            }
        }

        if (AppConstant.CHAR_EEG_CH3_SIGNAL == characteristic.uuid) {
            if (!mCh3!!.chEnabled) {
                mCh3!!.chEnabled = true
            }
            val mNewEEGdataBytes = characteristic.value
            val byteLength = mNewEEGdataBytes.size
            getDataRateBytes(byteLength)
            if (mEEGConnectedAllChannels) {
                mCh3!!.handleNewData(mNewEEGdataBytes, false)
                if (mCh3!!.packetCounter.toInt() == mPacketBuffer) {
                    addToGraphBuffer(mCh3!!, mGraphAdapterCh3, false)
                }
            }
        }

        if (AppConstant.CHAR_MPU_COMBINED == characteristic.uuid) {
            val bytes = characteristic.value
            getDataRateBytes2(bytes.size) //+=240
            mMPU!!.handleNewData(bytes, true)
//            addToGraphBufferMPU(mMPU!!)
            for (i in 0 until bytes!!.size / 12) {
                val ax = DataChannel.bytesToDoubleMPUAccel(bytes[12 * i], bytes[12 * i + 1])
                val ay = DataChannel.bytesToDoubleMPUAccel(bytes[12 * i + 2], bytes[12 * i + 3])
                mZAccValue = DataChannel.bytesToDoubleMPUAccel(bytes[12 * i + 4], bytes[12 * i + 5])
                val gx = DataChannel.bytesToDoubleMPUGyro(bytes[12 * i + 6], bytes[12 * i + 7])
                val gy = DataChannel.bytesToDoubleMPUGyro(bytes[12 * i + 8], bytes[12 * i + 9])
                val gz = DataChannel.bytesToDoubleMPUGyro(bytes[12 * i + 10], bytes[12 * i + 11])
                mSaveFileMPU?.exportDataDouble(ax, ay, mZAccValue, gx, gy, gz)
                // TODO : ADD TO BUFFER:
                mAccX?.addToBuffer(ax)
                mAccY?.addToBuffer(ay)
                mAccZ?.addToBuffer(mZAccValue)
                val zAccStr = "${mZAccValue.format(2)} g"
                runOnUiThread {
                    zAccelerationTextView.text = zAccStr
                }
            }
        }

        if (mCh1!!.chEnabled && mCh2!!.chEnabled && mCh3!!.chEnabled) {
            mNumberPackets++
            mEEGConnectedAllChannels = true
            mCh1!!.chEnabled = false
            mCh2!!.chEnabled = false
            mCh3!!.chEnabled = false
            if (mCh1!!.characteristicDataPacketBytes != null &&
                    mCh2!!.characteristicDataPacketBytes != null &&
                    mCh3!!.characteristicDataPacketBytes != null) {
                mPrimarySaveDataFile!!.writeToDisk(mCh1?.characteristicDataPacketBytes,
                        mCh2?.characteristicDataPacketBytes, mCh3?.characteristicDataPacketBytes)
            }
        }
    }

    private fun processClassifiedData(Y: Double) {
        //Shift backwards:
        System.arraycopy(yFitArray, 1, yFitArray, 0, 4)
        //Add to end;
        yFitArray[4] = Y
        //Analyze:
        Log.e(TAG, " YfitArray: " + Arrays.toString(yFitArray))
        val checkLastThreeMatches = lastThreeMatches(yFitArray)
        if (checkLastThreeMatches) {
            //Get value:
            Log.e(TAG, "Found fit: " + yFitArray[4].toString())
            val s = "[$Y]"
            runOnUiThread { mYfitTextView!!.text = s }
            // For Controlling Hand. Some commands have to be altered because the classes don't match the commands.
            // TODO: run BT Hand
            if (Y != 0.0) {
                val command: Int
                if (Y == 2.0) {
                    command = 1
                } else if (Y == 1.0) {
                    command = 2
                } else
                    command = Y.toInt()
                mConnectedThread.write(command)
            }
        } else {
            val b = (yFitArray[0] == 0.0 && yFitArray[1] == 0.0 && yFitArray[2] == 0.0
                    && yFitArray[3] == 0.0 && yFitArray[4] == 0.0)
            if (b) {
                val s = "[$Y]"
                runOnUiThread { mYfitTextView!!.text = s }
                // TODO: run BT Hand
                mConnectedThread.write(1)
            }
        }
    }

    private val mTrainEMGThread = Runnable {
        val dataConcat = Doubles.concat(mCh1!!.classificationBuffer, mCh2!!.classificationBuffer, mCh3!!.classificationBuffer, mEMGClassArray)
        // Assert Array Size.
        Log.e(TAG, "dataConcat.size: ${dataConcat.size}") // Should be 120k
        // Run feature extraction:
        mKNNParams = mNativeInterface.jTrainingRoutineKNN2(dataConcat)
        mClassifierReady = true
        Log.e(TAG, "mKNNPararms.size: ${mKNNParams.size}")
    }

    private val mClassifyEMGThread = Runnable {
        val arrayCh1 = DoubleArray(750)
        System.arraycopy(mCh1!!.classificationBuffer, 29249, arrayCh1, 0, 750)
        val arrayCh2 = DoubleArray(750)
        System.arraycopy(mCh2!!.classificationBuffer, 29249, arrayCh2, 0, 750)
        val arrayCh3 = DoubleArray(750)
        System.arraycopy(mCh3!!.classificationBuffer, 29249, arrayCh3, 0, 750)
        val concatAllChannels = Doubles.concat(arrayCh1, arrayCh2, arrayCh3)
        val yOutput = mNativeInterface.jClassifyUsingKNNv4(concatAllChannels, mKNNParams, 3.0)
        processClassifiedData(yOutput)
    }

//    private fun classifyEMG() {
//        if (mTFRunModel) {
//            val outputScores = FloatArray(4)
//            val end = mCh1!!.classificationBuffer.size - 1 //E.g. if size = 500, end = 499
//            val from = end - mTensorflowWindowSize
//            val ch1Input = Arrays.copyOfRange(mCh1!!.classificationBuffer, from, end)
//            val rescaledInput = mNativeInterface.jfiltRescale(ch1Input, mScaleFactor) // mTensorflowWindowSize
//            Log.i(TAG, "onCharacteristicChanged: TF_PRECALL_TIME, N#" + mNumberOfClassifierCalls.toString())
//            mTFInferenceInterface?.feed("keep_prob", floatArrayOf(1f))
//            mTFInferenceInterface?.feed(INPUT_DATA_FEED, rescaledInput, mTensorflowWindowSize.toLong())
//            mTFInferenceInterface?.run(arrayOf(OUTPUT_DATA_FEED))
//            mTFInferenceInterface?.fetch(OUTPUT_DATA_FEED, outputScores)
//            val yTF = DataChannel.getIndexOfLargest(outputScores)
//            Log.i(TAG, "CALL#" + mNumberOfClassifierCalls.toString() + ":\n" +
//                    "TF outputScores: " + Arrays.toString(outputScores))
//
//            System.arraycopy(mClassificationBuffer, 1, mClassificationBuffer, 0, 4)
//            mClassificationBuffer[4] = yTF
//            var output = 0
//            if (mMiniDrone != null) {
//                if (mMPUClass == 0) {
//                    if (allValuesSame(mClassificationBuffer)) {
//                        output = yTF
//                    }
//                } else {
//                    output = mMPUClass
//                }
//            } else {
//                output = yTF
//            }
//            sendDroneCommand(output)
//            executeWheelchairCommand(output)
//            Log.e(TAG, "Drone Class Output: $output")
//            mNumberOfClassifierCalls++
//            val s = Arrays.toString(mClassificationBuffer) + " - MPU: [$mMPUClass]"
//            val outputStr = "Output: [ $output ]"
//            runOnUiThread {
//                mYfitTextView!!.text = s
//                mEMGClassText?.text = outputStr
//            }
//        }
//    }

    private fun addToGraphBuffer(dataChannel: DataChannel, graphAdapter: GraphAdapter?, updateTrainingRoutine: Boolean) {
        if (mFilterData && dataChannel.totalDataPointsReceived > 1000) {
            graphAdapter!!.clearPlot()
            // TODO: Implement
        } else {
            if (dataChannel.dataBuffer != null) {
                if (mPrimarySaveDataFile!!.resolutionBits == 24) {
                    var i = 0
                    while (i < dataChannel.dataBuffer!!.size / 3) {
                        graphAdapter!!.addDataPointTimeDomain(DataChannel.bytesToDouble(dataChannel.dataBuffer!![3 * i],
                                dataChannel.dataBuffer!![3 * i + 1], dataChannel.dataBuffer!![3 * i + 2]),
                                dataChannel.totalDataPointsReceived - dataChannel.dataBuffer!!.size / 3 + i)
                        if (updateTrainingRoutine) {
                            for (j in 0 until graphAdapter.sampleRate / 250) {
                                updateTrainingRoutine(dataChannel.totalDataPointsReceived - dataChannel.dataBuffer!!.size / 3 + i + j)
                            }
                        }
                        i += graphAdapter.sampleRate / 250
                    }
                } else if (mPrimarySaveDataFile!!.resolutionBits == 16) {
                    var i = 0
                    while (i < dataChannel.dataBuffer!!.size / 2) {
                        graphAdapter!!.addDataPointTimeDomain(DataChannel.bytesToDouble(dataChannel.dataBuffer!![2 * i],
                                dataChannel.dataBuffer!![2 * i + 1]),
                                dataChannel.totalDataPointsReceived - dataChannel.dataBuffer!!.size / 2 + i)
                        if (updateTrainingRoutine) {
                            for (j in 0 until graphAdapter.sampleRate / 250) {
                                updateTrainingRoutine(dataChannel.totalDataPointsReceived - dataChannel.dataBuffer!!.size / 2 + i + j)
                            }
                        }
                        i += graphAdapter.sampleRate / 250
                    }
                }
            }
        }

        dataChannel.dataBuffer = null
        dataChannel.packetCounter = 0.toShort()
    }

    private fun updateTrainingRoutine(dataPoints: Int) {
        if (dataPoints % mSampleRate == 0 && mRunTrainingBool) {
            val second = dataPoints / mSampleRate
            val secondsBetweenAction = 10
            val countDownValue: Int
            if (second % secondsBetweenAction == 0) mMediaBeep.start()
            when {
                (second in 0 until 10) -> {
                    countDownValue = 10 - second
                    updateTrainingPromptColor(Color.GREEN)
                    mEMGClass = 0.0
                    updateTrainingPrompt("Relax Hand - Countdown to First Event: " + countDownValue + "s\n Next up: Close Hand")
                }
                (second in 10 until 20) -> {
                    countDownValue = 20 - second
                    mEMGClass = 1.0
                    updateTrainingPrompt("[1] Close Hand and Hold for " + countDownValue + "s\n Next up: Relax Hand", Color.RED)
                }
                (second in 20 until 30) -> {
                    countDownValue = 30 - second
                    mEMGClass = 0.0
                    updateTrainingPrompt("[0] Relax Hand and Remain for " + countDownValue + "s\n Next up: Close Pinky")
                }
                (second in 30 until 40) -> {
                    countDownValue = 40 - second
                    mEMGClass = 2.0
                    updateTrainingPrompt("[7] Close Pinky and Hold for " + countDownValue + "s\n Next up: Relax Hand", Color.RED)
                }
                (second in 40 until 50) -> {
                    countDownValue = 50 - second
                    mEMGClass = 0.0
                    updateTrainingPrompt("[0] Relax Hand and Remain for " + countDownValue + "s\n Next up: Close Ring")
                }
                (second in 50 until 60) -> {
                    countDownValue = 60 - second
                    mEMGClass = 3.0
                    updateTrainingPrompt("[6] Close Ring and Hold for " + countDownValue + "s\n Next up: Relax Hand", Color.RED)
                }
                (second in 60 until 70) -> {
                    countDownValue = 70 - second
                    mEMGClass = 0.0
                    updateTrainingPrompt("[0] Relax Hand and Remain for " + countDownValue + "s\n Next up: Close Middle")
                }
                (second in 70 until 80) -> {
                    countDownValue = 80 - second
                    mEMGClass = 4.0
                    updateTrainingPrompt("[5] Close Middle and Hold " + countDownValue + "s\n Next up: Relax Hand", Color.RED)
                }
                (second in 80 until 90) -> {
                    countDownValue = 90 - second
                    mEMGClass = 0.0
                    updateTrainingPrompt("[0] Relax Hand and Remain for " + countDownValue.toString() + "s\n Next up: Close Index")
                }
                (second in 90 until 100) -> {
                    countDownValue = 100 - second
                    mEMGClass = 5.0
                    updateTrainingPrompt("[4] Close Index and Hold " + countDownValue + "s\n Next up: Relax Hand", Color.RED)
                }
                (second in 100 until 110) -> {
                    countDownValue = 110 - second
                    mEMGClass = 0.0
                    updateTrainingPrompt("[0] Relax Hand and Remain for " + (countDownValue) + "s\n Next up: Close Thumb")
                }
                (second in 110 until 120) -> {
                    countDownValue = 120 - second
                    mEMGClass = 6.0
                    updateTrainingPrompt("[3] Close Thumb and Hold " + (countDownValue) + "s\n Next up: Relax Hand", Color.RED)
                }
                (second in 120 until 130) -> {
                    // TODO: Run training EXACTLY @ 120 seconds elapsed.
                    updateTrainingPrompt("Stop!")
                    updateTrainingPromptColor(Color.RED)
                    updateTrainingView(false)
                    mEMGClass = 0.0
                    mRunTrainingBool = false
                }
            }
        }

        if (dataPoints < 30000) {
            mEMGClassArray[dataPoints] = mEMGClass
            if (dataPoints == 29999) {
                // Do this once
                val trainThread = Thread(mTrainEMGThread)
                trainThread.start()
            }
        }
    }

    private fun updateTrainingPrompt(prompt: String, color: Int = Color.GREEN) {
        runOnUiThread {
            if (mRunTrainingBool) {
                mTrainingInstructions?.text = prompt
                mTrainingInstructions?.setTextColor(color)
            }
        }
    }

    private fun updateTrainingView(b: Boolean) {
        val visibility = if (b) View.VISIBLE else View.GONE
        runOnUiThread {
            mTrainingInstructions?.visibility = visibility
//            mEMGClassText?.visibility = visibility
        }
    }

    private fun updateTrainingPromptColor(color: Int) {
        runOnUiThread {
            if (mRunTrainingBool) {
                mTrainingInstructions?.setTextColor(color)
            }
        }
    }

    private fun getDataRateBytes2(bytes: Int) {
        val mCurrentTime = System.currentTimeMillis()
        points2 += bytes
        if (mCurrentTime > mLastTime2 + 3000) {
            val datarate2 = (points2 / 3).toDouble()
            points2 = 0
            mLastTime2 = mCurrentTime
            Log.e(" DataRate 2(MPU):", "$datarate2 Bytes/s")
        }
    }

    private fun getDataRateBytes(bytes: Int) {
        val mCurrentTime = System.currentTimeMillis()
        points += bytes
        if (mCurrentTime > mLastTime + 5000) {
            dataRate = (points / 5).toDouble()
            points = 0
            mLastTime = mCurrentTime
            Log.e(" DataRate:", "$dataRate Bytes/s")
            runOnUiThread {
                val s = "$dataRate Bytes/s"
                mDataRate!!.text = s
            }
        }
    }

    override fun onReadRemoteRssi(gatt: BluetoothGatt, rssi: Int, status: Int) {
        uiRssiUpdate(rssi)
    }

    override fun onConnectionStateChange(gatt: BluetoothGatt, status: Int, newState: Int) {
        when (newState) {
            BluetoothProfile.STATE_CONNECTED -> {
                mConnected = true
                runOnUiThread {
                    if (menu != null) {
                        menu!!.findItem(R.id.menu_connect).isVisible = false
                        menu!!.findItem(R.id.menu_disconnect).isVisible = true
                    }
                }
                Log.i(TAG, "Connected")
                updateConnectionState(getString(R.string.connected))
                invalidateOptionsMenu()
                runOnUiThread {
                    mDataRate!!.setTextColor(Color.BLACK)
                    mDataRate!!.setTypeface(null, Typeface.NORMAL)
                }
                //Start the service discovery:
                gatt.discoverServices()
                startMonitoringRssiValue()
            }
            BluetoothProfile.STATE_CONNECTING -> {
            }
            BluetoothProfile.STATE_DISCONNECTING -> {
            }
            BluetoothProfile.STATE_DISCONNECTED -> {
                mConnected = false
                runOnUiThread {
                    if (menu != null) {
                        menu!!.findItem(R.id.menu_connect).isVisible = true
                        menu!!.findItem(R.id.menu_disconnect).isVisible = false
                    }
                }
                Log.i(TAG, "Disconnected")
                runOnUiThread {
                    mDataRate!!.setTextColor(Color.RED)
                    mDataRate!!.setTypeface(null, Typeface.BOLD)
                    mDataRate!!.text = HZ
                }
                updateConnectionState(getString(R.string.disconnected))
                stopMonitoringRssiValue()
                invalidateOptionsMenu()
            }
            else -> {
            }
        }
    }

    private fun startMonitoringRssiValue() {
        readPeriodicallyRssiValue(true)
    }

    private fun stopMonitoringRssiValue() {
        readPeriodicallyRssiValue(false)
    }

    private fun readPeriodicallyRssiValue(repeat: Boolean) {
        mTimerEnabled = repeat
        // check if we should stop checking RSSI value
        if (!mConnected || !mTimerEnabled) {
            mTimerEnabled = false
            return
        }

        mTimerHandler.postDelayed(Runnable {
            if (!mConnected) {
                mTimerEnabled = false
                return@Runnable
            }
            // request RSSI value
            mBluetoothGattArray[0]!!.readRemoteRssi()
            // add call it once more in the future
            readPeriodicallyRssiValue(mTimerEnabled)
        }, RSSI_UPDATE_TIME_INTERVAL.toLong())
    }

    override fun onCharacteristicWrite(gatt: BluetoothGatt, characteristic: BluetoothGattCharacteristic, status: Int) {
        Log.i(TAG, "onCharacteristicWrite :: Status:: $status")
    }

    override fun onDescriptorWrite(gatt: BluetoothGatt, descriptor: BluetoothGattDescriptor, status: Int) {}

    override fun onDescriptorRead(gatt: BluetoothGatt, descriptor: BluetoothGattDescriptor, status: Int) {
        Log.i(TAG, "onDescriptorRead :: Status:: $status")
    }

    override fun onError(errorMessage: String) {
        Log.e(TAG, "Error:: $errorMessage")
    }

    private fun updateConnectionState(status: String) {
        runOnUiThread {
            if (status == getString(R.string.connected)) {
                Toast.makeText(applicationContext, "Device Connected!", Toast.LENGTH_SHORT).show()
            } else if (status == getString(R.string.disconnected)) {
                Toast.makeText(applicationContext, "Device Disconnected!", Toast.LENGTH_SHORT).show()
            }
        }
    }

    private fun updateBatteryStatus(integerValue: Int) {
        val status: String
        val convertedBatteryVoltage = integerValue.toDouble() / 4096.0 * 7.20
        //Because TPS63001 dies below 1.8V, we need to set up a linear fit between 1.8-4.2V
        //Anything over 4.2V = 100%
        val finalPercent: Double = when {
            125.0 / 3.0 * convertedBatteryVoltage - 75.0 > 100.0 -> 100.0
            125.0 / 3.0 * convertedBatteryVoltage - 75.0 < 0.0 -> 0.0
            else -> 125.0 / 3.0 * convertedBatteryVoltage - 75.0
        }
        Log.e(TAG, "Battery Integer Value: " + integerValue.toString())
        Log.e(TAG, "ConvertedBatteryVoltage: " + String.format(Locale.US, "%.5f", convertedBatteryVoltage) + "V : " + String.format(Locale.US, "%.3f", finalPercent) + "%")
        status = String.format(Locale.US, "%.1f", finalPercent) + "%"
        runOnUiThread {
            if (finalPercent <= batteryWarning) {
                mBatteryLevel!!.setTextColor(Color.RED)
                mBatteryLevel!!.setTypeface(null, Typeface.BOLD)
                Toast.makeText(applicationContext, "Charge Battery, Battery Low $status", Toast.LENGTH_SHORT).show()
            } else {
                mBatteryLevel!!.setTextColor(Color.GREEN)
                mBatteryLevel!!.setTypeface(null, Typeface.BOLD)
            }
            mBatteryLevel!!.text = status
        }
    }

    private fun uiRssiUpdate(rssi: Int) {
        runOnUiThread {
            val menuItem = menu!!.findItem(R.id.action_rssi)
            val statusActionItem = menu!!.findItem(R.id.action_status)
            val valueOfRSSI = rssi.toString() + " dB"
            menuItem.title = valueOfRSSI
            if (mConnected) {
                val newStatus = "Status: " + getString(R.string.connected)
                statusActionItem.title = newStatus
            } else {
                val newStatus = "Status: " + getString(R.string.disconnected)
                statusActionItem.title = newStatus
            }
        }
    }

    private val mMiniDroneListener = object : MiniDrone.Listener {
        override fun onDroneConnectionChanged(state: ARCONTROLLER_DEVICE_STATE_ENUM) {
            when (state) {
                ARCONTROLLER_DEVICE_STATE_ENUM.ARCONTROLLER_DEVICE_STATE_RUNNING -> mConnectionProgressDialog?.dismiss()

                ARCONTROLLER_DEVICE_STATE_ENUM.ARCONTROLLER_DEVICE_STATE_STOPPED -> {
                    // if the deviceController is stopped, go back to the previous activity
                    mConnectionProgressDialog?.dismiss()
                    finish()
                }

                else -> {
                }
            }
        }

        override fun onBatteryChargeChanged(batteryPercentage: Int) {
            val batteryFormatted = String.format(Locale.US, "Drone Battery: [%d%%]", batteryPercentage)
            mBatteryLevelDrone?.text = batteryFormatted
        }

        override fun onPilotingStateChanged(state: ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM) {
            when (state) {
                ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_LANDED -> {
                    val to = "Take off"
                    mTakeOffLandBt?.text = to
                    mTakeOffLandBt?.isEnabled = true
                }
                ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_FLYING, ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_ENUM.ARCOMMANDS_MINIDRONE_PILOTINGSTATE_FLYINGSTATECHANGED_STATE_HOVERING -> {
                    val ld = "Land"
                    mTakeOffLandBt?.text = ld
                    mTakeOffLandBt?.isEnabled = true
                }
                else -> mTakeOffLandBt?.isEnabled = false
            }
        }

        override fun onPictureTaken(error: ARCOMMANDS_MINIDRONE_MEDIARECORDEVENT_PICTUREEVENTCHANGED_ERROR_ENUM) {
            Log.i(TAG, "Picture has been taken")
        }

        override fun configureDecoder(codec: ARControllerCodec) {
            //            mVideoView.configureDecoder(codec)
        }

        override fun onFrameReceived(frame: ARFrame) {
            //            mVideoView.displayFrame(frame)
        }

        override fun onMatchingMediasFound(nbMedias: Int) {
            mDownloadProgressDialog?.dismiss()
            mNbMaxDownload = nbMedias
            mCurrentDownloadIndex = 1
            if (nbMedias > 0) {
                mDownloadProgressDialog = ProgressDialog(this@DeviceControlActivity, R.style.AppCompatAlertDialogStyle)
                mDownloadProgressDialog?.isIndeterminate = false
                mDownloadProgressDialog?.setProgressStyle(ProgressDialog.STYLE_HORIZONTAL)
                mDownloadProgressDialog?.setMessage("Downloading medias")
                mDownloadProgressDialog?.max = mNbMaxDownload * 100
                mDownloadProgressDialog?.secondaryProgress = mCurrentDownloadIndex * 100
                mDownloadProgressDialog?.progress = 0
                mDownloadProgressDialog?.setCancelable(false)
                mDownloadProgressDialog?.setButton(DialogInterface.BUTTON_NEGATIVE, "Cancel") { _, _ ->
                    mMiniDrone?.cancelGetLastFlightMedias()
                }
                mDownloadProgressDialog?.show()
            }
        }

        override fun onDownloadProgressed(mediaName: String, progress: Int) {
            mDownloadProgressDialog?.progress = (mCurrentDownloadIndex - 1) * 100 + progress
        }

        override fun onDownloadComplete(mediaName: String) {
            mCurrentDownloadIndex++
            mDownloadProgressDialog?.secondaryProgress = mCurrentDownloadIndex * 100

            if (mCurrentDownloadIndex > mNbMaxDownload) {
                mDownloadProgressDialog?.dismiss()
                mDownloadProgressDialog = null
            }
        }
    }

    private val mJSListener = object : JSDrone.Listener {
        override fun onDroneConnectionChanged(state: ARCONTROLLER_DEVICE_STATE_ENUM) {
            when (state) {
                ARCONTROLLER_DEVICE_STATE_ENUM.ARCONTROLLER_DEVICE_STATE_RUNNING -> mConnectionProgressDialog?.dismiss()

                ARCONTROLLER_DEVICE_STATE_ENUM.ARCONTROLLER_DEVICE_STATE_STOPPED -> {
                    // if the deviceController is stopped, go back to the previous activity
                    mConnectionProgressDialog?.dismiss()
                    finish()
                }

                else -> {
                }
            }
        }

        override fun onBatteryChargeChanged(batteryPercentage: Int) {
        }

        override fun onPictureTaken(error: ARCOMMANDS_JUMPINGSUMO_MEDIARECORDEVENT_PICTUREEVENTCHANGED_ERROR_ENUM) {
            Log.i(TAG, "Picture has been taken")
        }

        override fun configureDecoder(codec: ARControllerCodec) {}

        override fun onFrameReceived(frame: ARFrame) {
        }

        override fun onAudioStateReceived(inputEnabled: Boolean, outputEnabled: Boolean) {
        }

        override fun configureAudioDecoder(codec: ARControllerCodec) {
        }

        override fun onAudioFrameReceived(frame: ARFrame) {
        }

        override fun onMatchingMediasFound(nbMedias: Int) {
        }

        override fun onDownloadProgressed(mediaName: String, progress: Int) {
        }

        override fun onDownloadComplete(mediaName: String) {
        }
    }

    companion object {
        const val RECIEVE_MESSAGE = 1        // Status  for Handler
        const val HZ = "0 Hz"
        private val TAG = DeviceControlActivity::class.java.simpleName
        var mRedrawer: Redrawer? = null
        internal var mMPU: DataChannel? = null
        var mZAccValue = 1.0
        var mZAccMaxThreshold = 0.8
        var mZAccMinThreshold = 0.4
        internal var mAccX: DataChannel? = null
        internal var mAccY: DataChannel? = null
        internal var mAccZ: DataChannel? = null
        // Power Spectrum Graph Data:
        private var mSampleRate = 250
        var mEMGClass = 0.0
        //Data Channel Classes
        internal var mFilterData = false
        private var mPacketBuffer = 6
        //RSSI:
        private const val RSSI_UPDATE_TIME_INTERVAL = 2000
        //Save Data File
        private var mPrimarySaveDataFile: SaveDataFile? = null
        private var mSaveFileMPU: SaveDataFile? = null
        //Tensorflow CONSTANTS:
        // ADS1299 Register Configs
        val ADS1299_DEFAULT_BYTE_CONFIG = byteArrayOf(
                0x96.toByte(), 0xD0.toByte(), 0xEC.toByte(), 0x00.toByte(), //CONFIG1-3, LOFF
                0x40.toByte(), 0x40.toByte(), 0xE1.toByte(), 0xE1.toByte(), //CHSET 1-4
                0x00.toByte(), 0x00.toByte(), 0x00.toByte(), 0x00.toByte(), //CHSET 5-8
                0x0F.toByte(), 0x0F.toByte(), // BIAS_SENSP/N
                0x00.toByte(), 0x00.toByte(), 0x00.toByte(), 0x00.toByte(), 0x00.toByte(), // LOFF_P/N (IGNORE)
                0x0F.toByte(), 0x00.toByte(), 0x00.toByte(), 0x00.toByte()) //GPIO, MISC1 (0x20 for SRB1), MISC2, CONFIG4

        //Note for companion object: JNI call must include Companion in call: e.g. package_class_Companion_function(...).
        //TODO: Still does not work when I try to call from the companion object.

        // SPP UUID service
        val MY_UUID = UUID.fromString("00001101-0000-1000-8000-00805F9B34FB")!!

        // MAC-address of Bluetooth module (you must edit this line to change)
        const val address = "30:14:12:17:21:88"
        @JvmField
        var mHandler = Handler()

        init {
            System.loadLibrary("ssvep-lib")
        }
    }
}
